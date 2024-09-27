#![allow(unused)]
use std::collections::{HashMap, HashSet};
use cartesian::cartesian;

use itertools::Itertools;
use log::info;
use yui::bitseq::Bit;
use yui::{Ring, RingOps};
use yui_homology::{ChainComplexTrait, Grid1, XChainComplex, XModStr};
use yui_link::{Crossing, Edge, InvLink};

use crate::kh::v2::builder::TngComplexBuilder;
use crate::kh::v2::cob::LcCobTrait;
use crate::kh::v2::tng::{Tng, TngComp};
use crate::kh::v2::tng_complex::{TngComplex, TngKey};
use crate::kh::{KhChain, KhComplex, KhGen};
use crate::khi::{KhIChain, KhIComplex, KhIGen};

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    e_map: HashMap<Edge, Edge>,
    key_map: HashMap<TngKey, TngKey>,
    on_axis: HashSet<Crossing>,
    off_axis: HashSet<Crossing>,
    complex: TngComplex<R>,
    auto_validate: bool
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_kh_complex(l: &InvLink, h: &R, t: &R, reduced: bool, equivariant: bool) -> KhComplex<R> { 
        let mut b = Self::init(l, h, t, reduced);

        b.process_off_axis();
        b.process_on_axis_nonequiv();
        b.into_tng_complex().into_kh_complex(vec![])
    }

    pub fn build_khi_complex(l: &InvLink, h: &R, t: &R, reduced: bool, equivariant: bool) -> KhIComplex<R> { 
        let mut b = Self::init(l, h, t, reduced);

        b.process_off_axis();
        b.process_on_axis_equiv();
        b.finalize_equiv();;

        b.into_khi_complex()
    }

    pub fn init(l: &InvLink, h: &R, t: &R, reduced: bool) -> SymTngBuilder<R> { 
        assert!(l.link().is_knot(), "Only invertible knots are supported.");
        assert!(l.link().data().iter().all(|x| !x.is_resolved()));

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduced);
        let base_pt = if reduced { 
            // TODO check on-axis
            l.link().first_edge()
        } else { 
            None
        };

        let complex = TngComplex::init(h, t, deg_shift, base_pt);
        let e_map = l.link().edges().iter().map(|&e| (e, l.inv_e(e))).collect();
        let key_map = HashMap::new();
        let (on_axis, off_axis) = Self::separate_crossings(l);
        let auto_validate = false;

        SymTngBuilder { e_map, key_map, on_axis, off_axis, complex, auto_validate }
    }

    fn separate_crossings(l: &InvLink) -> (HashSet<Crossing>, HashSet<Crossing>) { 
        let mut on_axis = HashSet::new();
        let mut off_axis = HashSet::new();
        let mut take = false; // a flag to collect only half of the off-axis crossings.

        l.link().traverse_edges((0, 0), |i, j| { 
            let x = &l.link().data()[i];

            if i == l.inv_x(i) {
                take = !take;

                if !on_axis.contains(x) {
                    on_axis.insert(x.clone());
                }
            } else {
                let e = x.edge(j);
                if e == l.inv_e(e) { // on-axis edge
                    take = !take;
                }

                if take && !off_axis.contains(x) { 
                    off_axis.insert(x.clone());
                }
            }
        });

        assert_eq!(on_axis.len() + 2 * off_axis.len(), l.link().crossing_num());

        (on_axis, off_axis)
    }

    fn process_off_axis(&mut self) { 
        let (h, t) = self.complex.ht();

        let mut b = TngComplexBuilder::new(h, t, (0, 0), None); // half off-axis part
        b.set_crossings(self.off_axis.clone());
        b.process_all();

        let c1 = b.into_tng_complex();
        let c2 = c1.convert_edges(|e| self.inv_e(e));

        for ((k1, _), (k2, _)) in cartesian!(c1.iter_verts(), c2.iter_verts()) { 
            let k  = k1 + k2;
            let tk = k2 + k1;
            self.key_map.insert(k, tk);
        }

        self.complex.connect(c1);
        self.complex.connect(c2);

        if self.auto_validate { 
            self.validate_equiv();
        }
    }

    fn process_on_axis_nonequiv(&mut self) { 
        let complex = std::mem::take(&mut self.complex);
        
        let mut b = TngComplexBuilder::from(complex);
        b.set_crossings(self.on_axis.clone());
        b.process_all();

        self.complex = b.into_tng_complex();
    }

    fn process_on_axis_equiv(&mut self) { 
        let xs = self.on_axis.clone();

        for x in xs.iter() { 
            self.append_x(x);

            while let Some((k, r)) = self.find_good_equiv_loop(false) { 
                let updated = self.deloop_equiv(&k, r);
                for k in updated { 
                    if let Some((i, j)) = self.find_equiv_inv_edge(&k) { 
                        self.eliminate_equiv(&i, &j);
                    }
                }
            }
        }
    }

    fn append_x(&mut self, x: &Crossing) { 
        self.complex.append(x);

        if !x.is_resolved() { 
            // k <-> l  =>  k0 <-> l0,
            //              k1 <-> l1
            self.key_map = self.key_map.iter().flat_map(|(k, l)| { 
                let mut e0 = (*k, *l);
                let mut e1 = e0;
                e0.0.state.push(Bit::Bit0);
                e0.1.state.push(Bit::Bit0);
                e1.0.state.push(Bit::Bit1);
                e1.1.state.push(Bit::Bit1);
                [e0, e1]
            }).collect();
        }
    }

    fn find_loop(&self, allow_based: bool) -> Option<(TngKey, usize)> { 
        self.complex.iter_verts().find_map(|(k, v)| {
            v.tng().find_comp(|c| 
                c.is_circle() && 
                (allow_based || !self.complex.contains_base_pt(c))
            ).map(|r| (*k, r))
        })
    }

    fn find_good_equiv_loop(&self, allow_based: bool) -> Option<(TngKey, usize)> { 
        let find_in = |k: &TngKey, loc_tng: &Tng| { 
            if let Some(r_loc) = loc_tng.find_comp(|c| { 
                c.is_circle() && 
                (!self.is_sym_key(k) || self.is_sym_comp(c)) &&
                (allow_based || !self.complex.contains_base_pt(c))
            }) {
                let circ = loc_tng.comp(r_loc);
                let v = self.complex.vertex(k);
                let r = v.tng().find_comp(|c| c == circ).unwrap();
                Some((*k, r))
            } else {
                None
            }
        };
        
        self.complex.iter_verts().find_map(|(k, v)|
            v.out_edges().filter(|l| 
                self.is_equiv_edge(k, l)
            ).find_map(|l|
                self.complex.edge(k, l).gens().find_map(|cob| {
                    cob.comps().find_map(|c| 
                        if c.is_cap() || c.is_merge() { 
                            find_in(k, c.src())
                        } else if c.is_cup() || c.is_split() {
                            find_in(l, c.tgt())
                        } else { 
                            None
                        }
                    )
                })
            )
        )
    }

    fn deloop_equiv(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> { 
        let c = self.complex.vertex(&k).tng().comp(r);

        let updated = if self.is_sym_key(&k) { 
            if self.is_sym_comp(c) {
                // symmetric loop on symmetric key
                self.deloop_sym(&k, r)
            } else {
                // asymmetric loop on symmetric key
                self.deloop_asym(&k, r)
            }
        } else { 
            // (symmetric or asymmetric) loop on asymmetric key
            self.deloop_parallel(&k, r)
        };

        if self.auto_validate { 
            self.validate_equiv();
        }

        updated
    }

    fn deloop_sym(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex.vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(self.is_sym_comp(c));

        info!("deloop sym: {c} in {}", self.complex.vertex(&k));

        let updated = self.complex.deloop(k, r);

        self.remove_key_pair(k);

        for &k_new in updated.iter() { 
            self.add_key_pair(k_new, k_new);
        }

        updated
    }

    #[allow(non_snake_case)]
    fn deloop_asym(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex.vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(!self.is_sym_comp(c));

        // TODO support deloop on based loop
        assert!(!self.complex.contains_base_pt(c));

        //          ⚪︎1 | ⚪︎1
        //  ⚪︎1 | ⚪︎X  <-->  ⚪︎X | ⚪︎1
        //          ⚪︎X | ⚪︎X

        let tc = c.convert_edges(|e| self.inv_e(e));

        info!("deloop asym: {c} <--> {tc} in {}", self.complex.vertex(&k));

        let ks = self.complex.deloop(k, r);

        let (k_X, k_1) = (ks[0], ks[1]);
        let (k_XX, k_X1) = { 
            let tr = self.complex.vertex(&k_X).tng().index_of(&tc).unwrap();
            let tks = self.complex.deloop(&k_X, tr);
            (tks[0], tks[1])
        };
        let (k_1X, k_11) = { 
            let tr = self.complex.vertex(&k_1).tng().index_of(&tc).unwrap();
            let tks = self.complex.deloop(&k_1, tr);
            (tks[0], tks[1])
        };

        self.remove_key_pair(k);

        self.add_key_pair(k_XX, k_XX);
        self.add_key_pair(k_X1, k_1X);
        self.add_key_pair(k_11, k_11);

        vec![k_XX, k_X1, k_11]
    }

    #[allow(non_snake_case)]
    fn deloop_parallel(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex.vertex(k).tng().comp(r);

        assert!(!self.is_sym_key(k));

        //  ⚪︎1 | ..  <-->  .. | ⚪︎1
        //  ⚪︎X | ..  <-->  .. | ⚪︎X

        let tk = *self.inv_key(k);
        let tc = c.convert_edges(|e| self.inv_e(e));
        let tr = self.complex.vertex(&tk).tng().index_of(&tc).unwrap();

        info!("deloop parallel:");
        info!("  {c} in {}", self.complex.vertex(&k));
        info!("  {tc} in {}", self.complex.vertex(&tk));
        
        let ks = self.complex.deloop(k, r);
        let tks = self.complex.deloop(&tk, tr);

        self.remove_key_pair(k);

        for (&k_new, &tk_new) in Iterator::zip(ks.iter(), tks.iter()) { 
            self.add_key_pair(k_new, tk_new);
        }

        ks
    }

    pub fn find_equiv_inv_edge(&self, k: &TngKey) -> Option<(TngKey, TngKey)> { 
        let v = self.complex.vertex(k);
        let sym = self.is_sym_key(k);
        let edges = Iterator::chain(
            v.in_edges().filter_map(|i| 
                (self.is_sym_key(i) == sym).then_some((i, k))
            ),
            v.out_edges().filter_map(|i| 
                (self.is_sym_key(i) == sym).then_some((k, i))
            )
        );

        edges.filter_map(|(i, j)| {
            self.is_equiv_inv_edge(i, j).then(|| {
                let s = self.complex.elim_weight(i, j);
                (s, i, j)
            })
        })
        .min_by_key(|(s, _, _)| *s)
        .map(|(_, i, j)| (*i, *j))
    }

    fn is_equiv_inv_edge(&self, i: &TngKey, j: &TngKey) -> bool { 
        let f = self.complex.edge(i, j);
        f.is_invertible() && self.is_equiv_edge(i, j)
    }

    fn is_equiv_edge(&self, i: &TngKey, j: &TngKey) -> bool { 
        if self.is_sym_key(i) && self.is_sym_key(j) { 
            true
        } else if !self.is_sym_key(i) && !self.is_sym_key(j) { 
            //  i - - -> j 
            //    \   /   
            //      /     : not allowed
            //    /   \   
            // ti - - -> tj
            let ti = self.inv_key(i);
            let tj = self.inv_key(j);

            !self.complex.keys_into(j).contains(ti) && 
            !self.complex.keys_into(tj).contains(i)
        } else { 
            false
        }
    }

    pub fn eliminate_equiv(&mut self, i: &TngKey, j: &TngKey) {
        assert_eq!(self.is_sym_key(i), self.is_sym_key(j));

        if self.is_sym_key(i) { 
            info!("eliminate sym: {i} -> {j}: {}", self.complex.edge(i, j));
            self.complex.eliminate(i, j);
        } else { 
            let ti = *self.inv_key(i);
            let tj = *self.inv_key(j);

            info!("eliminate parallel:");
            info!("  {i} -> {j}: {}", self.complex.edge(i, j));
            info!("  {ti} -> {tj}: {}", self.complex.edge(&ti, &tj));

            assert!(self.complex.has_edge(&i, &j));
            assert!(self.complex.has_edge(&ti, &tj));

            self.complex.eliminate(i, j);

            assert!(self.complex.has_edge(&ti, &tj));
            self.complex.eliminate(&ti, &tj);
        }

        self.remove_key_pair(i);
        self.remove_key_pair(j);

        if self.auto_validate { 
            self.validate_equiv();
        }
    }

    fn finalize_equiv(&mut self) {
        info!("finalize");

        for b in [false, true] { 
            while let Some((k, r)) = self.find_loop(b) { 
                self.deloop_equiv(&k, r);
            }
        }    
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }

    fn into_khi_complex(mut self) -> KhIComplex<R> {
        assert!(self.complex.is_finalizable());
        
        let n = self.complex.dim();
        let deg_shift = self.complex.deg_shift();
        let i0 = deg_shift.0;
        let i1 = i0 + (n as isize) + 1;

        let key_map = std::mem::take(&mut self.key_map);
        let inv_key = move |k: &TngKey| -> TngKey { 
            key_map[k]
        };

        let c = self.complex.into_raw_complex();

        let summands = Grid1::generate(i0..=i1, |i| { 
            let b_gens = c[i].gens().iter().map(|x| KhIGen::B(*x));
            let q_gens = c[i - 1].gens().iter().map(|x| KhIGen::Q(*x));
            XModStr::free(b_gens.chain(q_gens))
        });

        let d = move |i: isize, x: &KhIGen| -> KhIChain<R> { 
            match x { 
                KhIGen::B(x) => {
                    let z = KhChain::from(*x);
                    let dx = c.d(i, &z).map_gens(|y| KhIGen::B(*y));
                    let qx = KhIChain::from(KhIGen::Q(*x));
                    let qtx = {
                        let k = TngKey::from(x);
                        let tk = inv_key(&k);
                        let tx = tk.as_gen(deg_shift);
                        KhIChain::from(KhIGen::Q(tx))
                    };
                    dx + qx + qtx
                },
                KhIGen::Q(x) => {
                    let z = KhChain::from(*x);
                    c.d(i, &z).map_gens(|y| KhIGen::Q(*y))
                }
            }
        };

        let inner = XChainComplex::new(summands, 1, move |i, z| { 
            z.apply(|x| d(i, x))
        });

        KhIComplex::new_impl(inner, vec![], deg_shift)
    }

    fn inv_e(&self, e: Edge) -> Edge { 
        self.e_map[&e]
    }

    fn inv_key(&self, k: &TngKey) -> &TngKey { 
        &self.key_map[&k]
    }

    fn add_key_pair(&mut self, k: TngKey, tk: TngKey) { 
        self.key_map.insert(k, tk);
        if k != tk { 
            self.key_map.insert(tk, k);
        }
    }

    fn remove_key_pair(&mut self, k: &TngKey) { 
        let tk = self.key_map.remove(k).unwrap();
        if k != &tk { 
            self.key_map.remove(&tk);
        }
    }

    fn is_sym_key(&self, k: &TngKey) -> bool { 
        self.inv_key(k) == k
    }

    fn is_sym_comp(&self, c: &TngComp) -> bool { 
        &c.convert_edges(|e| self.inv_e(e)) == c
    }

    #[allow(unused)]
    fn print_keys(&self) {
        let mut done = HashSet::new();
        for k in self.key_map.keys().sorted() { 
            if done.contains(&k) { continue }

            let tk = self.inv_key(k);
            if k == tk {
                println!("{}", self.complex.vertex(k));
            } else { 
                println!("{} ↔ {}", self.complex.vertex(k), self.complex.vertex(tk));
            }

            done.insert(k);
            done.insert(tk);
        }
        println!();
    }

    fn validate_equiv(&self) {
        for (k, _) in self.complex.iter_verts() { 
            assert!(self.key_map.contains_key(k), "no inv-key for {k}");
            let tk = self.inv_key(k);

            for l in self.complex.keys_out_from(k) { 
                assert!(self.key_map.contains_key(l), "no inv-key for {l}");
                let tl = self.inv_key(l);

                assert!(self.complex.has_edge(tk, tl));

                let f = self.complex.edge(k, l);
                let tf = self.complex.edge(tk, tl);

                assert_eq!(&f.convert_edges(|e| self.inv_e(e)), tf);
            }
        }
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests { 
    use super::*;
    use num_traits::Zero;

    use yui::FF2;
    use yui::poly::HPoly;
    use yui_homology::{ChainComplexCommon, DisplaySeq, DisplayTable, RModStr};

    fn init_logger(l: log::LevelFilter) {
        use simplelog::*;
        TermLogger::init(
            l,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto
        ).unwrap()
    }

    #[test]
    fn test_kh_3_1() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());
        let c = SymTngBuilder::build_kh_complex(&l, &h, &t, false, true);
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn test() { 
        type R = HPoly<'H', FF2>;
    
        init_logger(log::LevelFilter::Debug);

        let l = InvLink::load("5_1").unwrap();
        let (h, t) = (R::variable(), R::zero());
        let _c = SymTngBuilder::build_khi_complex(&l, &h, &t, false, true);
    }
}