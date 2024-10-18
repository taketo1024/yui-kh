use std::collections::{HashMap, HashSet};
use std::ops::RangeInclusive;
use cartesian::cartesian;

use itertools::Itertools;
use log::info;
use yui::bitseq::Bit;
use yui::{RangeExt, Ring, RingOps};
use yui_link::{Crossing, Edge, InvLink};

use crate::kh::{KhComplex, KhGen};
use crate::khi::KhIComplex;
use crate::kh::internal::v2::builder::TngComplexBuilder;
use crate::kh::internal::v2::cob::LcCobTrait;
use crate::kh::internal::v2::tng::TngComp;
use crate::kh::internal::v2::tng_complex::{TngComplex, TngKey};

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: TngComplexBuilder<R>,
    e_map: HashMap<Edge, Edge>,
    key_map: HashMap<TngKey, TngKey>,
    on_axis: HashSet<Crossing>,
    off_axis: HashSet<Crossing>,
    auto_validate: bool
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_kh_complex(l: &InvLink, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> KhComplex<R> { 
        let mut b = Self::new(l, h, t, reduced);
        if let Some(h_range) = h_range { 
            b.inner.set_h_range(h_range);
        }
        b.process_all();
        b.into_kh_complex()
    }

    pub fn build_khi_complex(l: &InvLink, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> KhIComplex<R> { 
        let mut b = Self::new(l, h, t, reduced);
        if let Some(h_range) = h_range { 
            b.inner.set_h_range(h_range.mv(-1, 0)); // must extend left
        }
        b.process_all();
        b.into_khi_complex()
    }

    pub fn new(l: &InvLink, h: &R, t: &R, reduced: bool) -> SymTngBuilder<R> { 
        assert!(l.link().is_knot(), "Only invertible knots are supported.");
        assert!(l.link().data().iter().all(|x| !x.is_resolved()));
        assert!(!reduced || l.base_pt().is_some());

        let base_pt = if reduced { l.base_pt() } else { None };
        let mut inner = TngComplexBuilder::new(l.link(), h, t, base_pt);
        inner.set_crossings(vec![]); // clear crossings

        let e_map = l.link().edges().iter().map(|&e| (e, l.inv_e(e))).collect();
        let key_map = HashMap::new();
        let (on_axis, off_axis) = Self::separate_crossings(l);

        let auto_validate = false;

        SymTngBuilder { inner, e_map, key_map, on_axis, off_axis, auto_validate }
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

    pub fn process_all(&mut self) { 
        self.process_off_axis();
        self.process_on_axis();
        self.finalize();
    }

    fn process_off_axis(&mut self) { 
        info!("process off-axis: [{}]", self.off_axis.iter().join(", "));

        // process half
        let off_axis = std::mem::take(&mut self.off_axis);
        self.inner.set_crossings(off_axis);
        self.inner.process_all();

        // preserve old keys
        let keys = self.complex().keys().cloned().collect_vec();

        info!("merge two sides...");
        info!("  current: {}", self.inner.stat());

        // create other half and connect
        let mut tc = self.complex().convert_edges(|e| self.inv_e(e));
        tc.set_deg_shift((0, 0));
        
        self.inner.complex_mut().connect(tc);

        info!("     done: {}", self.inner.stat());

        // create other half for each element
        let inv_e = self.e_map.clone();
        self.inner.elements_mut().iter_mut().for_each(|e|
            e.modify(|k, cob| (
                k + k, 
                cob.map_cob(|c| {
                    let tc = c.convert_edges(|e| inv_e[&e]);
                    c.connect(tc)
                })
            ))
        );

        // update keys
        for (k1, k2) in cartesian!(keys.iter(), keys.iter()) { 
            let k  = k1 + k2;
            let tk = k2 + k1;
            self.key_map.insert(k, tk);
        }

        // drop if possible
        self.inner.drop_vertices();

        if self.auto_validate { 
            self.validate_equiv();
        }
    }

    fn process_on_axis(&mut self) { 
        info!("process on-axis: [{}]", self.on_axis.iter().join(", "));
        
        let on_axis = std::mem::take(&mut self.on_axis);
        self.inner.set_crossings(on_axis);
        self.inner.auto_deloop = false;
        self.inner.auto_elim = false;

        while let Some(x) = self.inner.choose_next() { 
            self.inner.process(&x);
            self.append_x(&x);

            while let Some((k, r)) = self.inner.find_good_loop(false) { 
                self.deloop_equiv(&k, r);
            }

            while let Some((k, r)) = self.inner.find_loop(false) { 
                self.deloop_equiv(&k, r);
            }
        }
    }

    fn append_x(&mut self, x: &Crossing) { 
        if x.is_resolved() { return } 

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

    fn deloop_equiv(&mut self, k: &TngKey, r: usize) { 
        let c = self.complex().vertex(&k).tng().comp(r);

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

        for k in updated { 
            if let Some((i, j)) = self.find_equiv_inv_edge(&k) { 
                self.eliminate_equiv(&i, &j);
            }
        }
    }

    fn deloop_sym(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex().vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(self.is_sym_comp(c));

        let updated = self.inner.deloop(k, r);

        self.remove_key_pair(k);

        for &k_new in updated.iter() { 
            self.add_key_pair(k_new, k_new);
        }

        updated
    }

    #[allow(non_snake_case)]
    fn deloop_asym(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex().vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(!self.is_sym_comp(c));
        assert!(!self.complex().contains_base_pt(c));

        //          ⚪︎1 | ⚪︎1
        //  ⚪︎1 | ⚪︎X  <-->  ⚪︎X | ⚪︎1
        //          ⚪︎X | ⚪︎X

        let tc = c.convert_edges(|e| self.inv_e(e));

        let ks = self.inner.deloop(k, r);

        let (k_X, k_1) = (ks[0], ks[1]);
        let (k_XX, k_X1) = { 
            let tr = self.complex().vertex(&k_X).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_X, tr);
            (tks[0], tks[1])
        };
        let (k_1X, k_11) = { 
            let tr = self.complex().vertex(&k_1).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_1, tr);
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
        let c = self.complex().vertex(k).tng().comp(r);

        assert!(!self.is_sym_key(k));

        //  ⚪︎1 | ..  <-->  .. | ⚪︎1
        //  ⚪︎X | ..  <-->  .. | ⚪︎X

        let tk = *self.inv_key(k);
        let tc = c.convert_edges(|e| self.inv_e(e));
        let tr = self.complex().vertex(&tk).tng().index_of(&tc).unwrap();

        let ks = self.inner.deloop(k, r);
        let tks = self.inner.deloop(&tk, tr);

        self.remove_key_pair(k);

        for (&k_new, &tk_new) in Iterator::zip(ks.iter(), tks.iter()) { 
            self.add_key_pair(k_new, tk_new);
        }

        ks
    }

    pub fn find_equiv_inv_edge(&self, k: &TngKey) -> Option<(TngKey, TngKey)> { 
        let v = self.complex().vertex(k);
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
                let s = self.complex().elim_weight(i, j);
                (s, i, j)
            })
        })
        .min_by_key(|(s, _, _)| *s)
        .map(|(_, i, j)| (*i, *j))
    }

    fn is_equiv_inv_edge(&self, i: &TngKey, j: &TngKey) -> bool { 
        let f = self.complex().edge(i, j);
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

            !self.complex().keys_into(j).contains(ti) && 
            !self.complex().keys_into(tj).contains(i)
        } else { 
            false
        }
    }

    pub fn eliminate_equiv(&mut self, i: &TngKey, j: &TngKey) {
        assert_eq!(self.is_sym_key(i), self.is_sym_key(j));
        assert!(self.complex().has_edge(&i, &j));

        if self.is_sym_key(i) { 
            self.inner.eliminate(i, j);
        } else { 
            let ti = *self.inv_key(i);
            let tj = *self.inv_key(j);

            assert!(self.complex().has_edge(&ti, &tj));

            self.inner.eliminate(i, j);
            self.inner.eliminate(&ti, &tj);
        }

        self.remove_key_pair(i);
        self.remove_key_pair(j);

        if self.auto_validate { 
            self.validate_equiv();
        }
    }

    fn finalize(&mut self) {
        info!("finalize");

        for b in [false, true] { 
            while let Some((k, r)) = self.inner.find_loop(b) { 
                self.deloop_equiv(&k, r);
            }
        }
    }

    pub fn complex(&self) -> &TngComplex<R> { 
        self.inner.complex()
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.inner.into_tng_complex()
    }

    pub fn into_kh_complex(self) -> KhComplex<R> {
        self.inner.into_kh_complex()
    }

    fn into_khi_complex(mut self) -> KhIComplex<R> {
        assert!(self.complex().is_finalizable());
        
        let deg_shift = self.complex().deg_shift();
        let h_range = self.inner.h_range();
        let key_map = std::mem::take(&mut self.key_map);

        let map = move |x: &KhGen| -> KhGen { 
            let k = TngKey::from(x);
            let tk = key_map[&k];
            let tx = tk.as_gen(deg_shift);
            tx
        };

        let c = self.into_kh_complex();
        let ci = KhIComplex::from_kh_complex(c, map);

        if let Some(h_range) = h_range { 
            ci.truncated(h_range.mv(1, 0))
        } else { 
            ci
        }
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
                println!("{}", self.complex().vertex(k));
            } else { 
                println!("{} ↔ {}", self.complex().vertex(k), self.complex().vertex(tk));
            }

            done.insert(k);
            done.insert(tk);
        }
        println!();
    }

    fn validate_equiv(&self) {
        for (k, _) in self.complex().iter_verts() { 
            assert!(self.key_map.contains_key(k), "no inv-key for {k}");
            let tk = self.inv_key(k);

            for l in self.complex().keys_out_from(k) { 
                assert!(self.key_map.contains_key(l), "no inv-key for {l}");
                let tl = self.inv_key(l);

                assert!(self.complex().has_edge(tk, tl));

                let f = self.complex().edge(k, l);
                let tf = self.complex().edge(tk, tl);

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

    #[test]
    fn test_kh_3_1() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());

        let c = SymTngBuilder::build_kh_complex(&l, &h, &t, false, None);
        c.check_d_all();;

        let h = c.inner().homology(false);

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn test_khi_3_1() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());

        let c = SymTngBuilder::build_khi_complex(&l, &h, &t, false, None);
        c.check_d_all();;

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 4);
        assert_eq!(h[4].rank(), 2);
    }

    #[test]
    fn test_kh_6_3_partial() { 
        let l = InvLink::load("6_3").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());
        let range = -1..=2;

        let c = SymTngBuilder::build_kh_complex(&l, &h, &t, false, Some(range.clone()));
        assert_eq!(c.h_range(), range);
        c.check_d_all();

        let h = c.homology();
        assert_eq!(h[0].rank(), 6);
        assert_eq!(h[1].rank(), 4);
    }
    
    #[test]
    fn test_khi_6_3_partial() { 
        let l = InvLink::load("6_3").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());
        let range = -1..=2;

        let c = SymTngBuilder::build_khi_complex(&l, &h, &t, false, Some(range.clone()));
        assert_eq!(c.h_range(), range);
        c.check_d_all();

        let h = c.homology();
        assert_eq!(h[0].rank(), 10);
        assert_eq!(h[1].rank(), 10);
    }
}