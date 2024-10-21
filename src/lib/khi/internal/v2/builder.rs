use std::collections::{HashMap, HashSet};
use std::ops::RangeInclusive;
use cartesian::cartesian;

use itertools::Itertools;
use log::info;
use yui::bitseq::Bit;
use yui::{KeyedUnionFind, RangeExt, Ring, RingOps};
use yui_link::{Crossing, Edge, InvLink};

use crate::kh::{KhComplex, KhGen};
use crate::khi::KhIComplex;
use crate::kh::internal::v2::builder::{BuildElem, TngComplexBuilder};
use crate::kh::internal::v2::cob::LcCobTrait;
use crate::kh::internal::v2::tng::TngComp;
use crate::kh::internal::v2::tng_complex::{TngComplex, TngKey};

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: TngComplexBuilder<R>,
    x_map: HashMap<Crossing, Crossing>,
    e_map: HashMap<Edge, Edge>,
    key_map: HashMap<TngKey, TngKey>,
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_kh_complex(l: &InvLink, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> KhComplex<R> { 
        let mut b = Self::new(l, h, t, reduced);
        if let Some(h_range) = h_range { 
            b.inner.set_h_range(h_range);
        }
        b.process_all();
        b.finalize();
        b.into_kh_complex()
    }

    pub fn build_khi_complex(l: &InvLink, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> KhIComplex<R> { 
        let mut b = Self::new(l, h, t, reduced);
        if let Some(h_range) = h_range { 
            b.inner.set_h_range(h_range.mv(-1, 0)); // must extend left
        }
        b.process_all();
        b.finalize();
        b.into_khi_complex()
    }

    pub fn new(l: &InvLink, h: &R, t: &R, reduced: bool) -> SymTngBuilder<R> { 
        assert!(l.link().data().iter().all(|x| !x.is_resolved()));
        assert!(!reduced || l.base_pt().is_some());

        let base_pt = if reduced { l.base_pt() } else { None };

        let mut inner = TngComplexBuilder::new(l.link(), h, t, base_pt);
        inner.auto_deloop = false;
        inner.auto_elim = false;

        let x_map = l.link().data().iter().map(|x| (x.clone(), l.inv_x(x).clone())).collect();
        let e_map = l.link().edges().iter().map(|&e| (e, l.inv_e(e))).collect();
        let key_map = HashMap::from_iter([(TngKey::init(), TngKey::init())]);

        SymTngBuilder { inner, x_map, e_map, key_map }
    }

    pub fn complex(&self) -> &TngComplex<R> { 
        self.inner.complex()
    }

    pub fn set_crossings<I>(&mut self, crossings: I) 
    where I: IntoIterator<Item = Crossing>{ 
        self.inner.set_crossings(crossings);
    }

    pub fn set_elements<I>(&mut self, elements: I) 
    where I: IntoIterator<Item = BuildElem<R>>{ 
        self.inner.set_elements(elements);
    }
    
    pub fn process_all(&mut self) { 
        self.process_off_axis();
        self.process_on_axis();
    }

    fn process_off_axis(&mut self) { 
        assert_eq!(self.complex().dim(), 0, "must start from init state.");

        info!("process off-axis");

        let off_axis = self.extract_off_axis_crossings(true);
        let elements = self.inner.take_elements();
        let b = self.build_from_half(off_axis, elements);

        self.merge(b);

        // drop if possible
        self.inner.drop_vertices();
    }

    fn extract_off_axis_crossings(&mut self, take_half: bool) -> Vec<Crossing> { 
        let crossings = self.inner.take_crossings();
        let (on_axis, off_axis) = crossings.into_iter().partition::<Vec<_>, _>(|x|
            self.inv_x(x) == x
        );
        self.inner.set_crossings(on_axis);

        if !take_half { 
            return off_axis
        }

        let is_adj = |x: &Crossing, y: &Crossing| {
            x.edges().iter().filter(|&&e| 
                self.inv_e(e) != e // no axis-crossing edge
            ).any(|e| 
                y.edges().contains(e)
            )
        };

        let mut u = KeyedUnionFind::from_iter(off_axis.iter());

        for (i, x) in off_axis.iter().enumerate() { 
            for j in 0 .. i { 
                let y = &off_axis[j];
                if is_adj(x, y) { 
                    u.union(&x, &y);
                }
            }
        }

        let retain = u.group().into_iter().fold(HashSet::<&Crossing>::new(), |mut res, next| { 
            if let Some(x) = next.iter().next() { 
                let tx = self.inv_x(x);
                if !res.contains(&tx) { 
                    res.extend(next.into_iter());
                }
            }
            res
        });

        let flags = off_axis.iter().map(|x| retain.contains(x)).collect_vec();

        off_axis.into_iter().enumerate().filter(|(i, _)| 
            flags[*i]
        ).map(|(_, x)| x).collect()
    }

    fn build_from_half(&self, crossings: Vec<Crossing>, elements: Vec<BuildElem<R>>) -> Self { 
        let (h, t) = self.complex().ht();
        let mut b = TngComplexBuilder::init(h, t, (0, 0), None);

        b.set_crossings(crossings);
        b.set_elements(elements);
        b.process_all();

        // take results
        let keys = b.complex().keys().cloned().collect_vec();
        let elements = b.take_elements();
        let mut c = b.into_tng_complex();

        // merge complexes
        info!("merge two sides...");
        info!("  current: {}", c.stat());

        let tc = c.convert_edges(|e| self.inv_e(e));
        c.connect(tc);

        info!("     done: {}", c.stat());

        // build elements
        let elements = elements.into_iter().map(|mut e| { 
            e.modify(|k, c| {
                let kk = k + k;
                let cc = c.into_map(|mut c, r| { 
                    let r = &r * &r;
                    let tc = c.convert_edges(|e| self.inv_e(e));
                    c.connect(tc);
                    (c, r)
                });
                (kk, cc)
            });
            e
        }).collect_vec();

        // build keys
        let key_map = cartesian!(keys.iter(), keys.iter()).map(|(k1, k2)| { 
            let k  = k1 + k2;
            let tk = k2 + k1;
            (k, tk)
        }).collect::<HashMap<_, _>>();

        let mut inner = TngComplexBuilder::from(c);
        inner.set_elements(elements);

        Self {
            inner,
            x_map: self.x_map.clone(),
            e_map: self.e_map.clone(),
            key_map,
        }
    }

    fn process_on_axis(&mut self) { 
        info!("process on-axis");
        
        while let Some(x) = self.inner.choose_next() { 
            self.append_x(&x);
            self.deloop_all(false);
        }
    }

    fn append_x(&mut self, x: &Crossing) { 
        self.inner.process(&x);
        
        if x.is_resolved() { return } 

        // k <-> l  =>  k0 <-> l0,
        //              k1 <-> l1
        self.key_map = self.key_map.iter().flat_map(|(k, l)| { 
            [Bit::Bit0, Bit::Bit1].map(|b| { 
                let mut k = *k;
                let mut l = *l;
                k.state.push(b);
                l.state.push(b);
                (k, l)
            })
        }).collect();
    }

    pub fn deloop_all(&mut self, allow_based: bool) { 
        while let Some((k, r)) = self.inner.find_good_loop(allow_based) { 
            self.deloop_equiv(&k, r);
        }
        while let Some((k, r)) = self.inner.find_loop(allow_based) { 
            self.deloop_equiv(&k, r);
        }
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
    }

    pub fn finalize(&mut self) {
        info!("finalize");

        self.deloop_all(false);
        self.deloop_all(true);

        assert!(self.complex().is_completely_delooped());
    }

    pub fn merge(&mut self, mut other: Self) { 
        let key_map = std::mem::take(&mut other.key_map);
        let elements = other.inner.take_elements();
        let c = other.into_tng_complex();

        self.inner.complex_mut().connect(c);
        self.key_map = self.key_map.iter().flat_map(|(k1, l1)|
            key_map.iter().map(move |(k2, l2)|
                (k1 + k2, l1 + l2)
            )
        ).collect();

        // TODO merge elements
        self.inner.set_elements(elements);
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.inner.into_tng_complex()
    }

    pub fn into_kh_complex(self) -> KhComplex<R> {
        self.inner.into_kh_complex()
    }

    pub fn into_khi_complex(mut self) -> KhIComplex<R> {
        assert!(self.complex().is_completely_delooped());
        
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

    fn inv_x(&self, x: &Crossing) -> &Crossing { 
        &self.x_map[x]
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

    #[allow(unused)]
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

    pub fn stat(&self) -> String { 
        self.inner.stat()
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests { 
    use crate::khi::KhIHomology;

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

    #[test]
    fn test_merge() { 
        let l = InvLink::sinv_knot_from_code([
            [6,9,7,10],[8,1,9,2],[14,7,1,8], // upper
            [3,13,4,12],[10,5,11,6],[11,3,12,2],[13,5,14,4], // lower
        ]); // 6_3

        let (h, t) = (FF2::zero(), FF2::zero());
        let mut b1 = SymTngBuilder::new(&l, &h, &t, false);
        let mut b2 = SymTngBuilder::new(&l, &h, &t, false);

        b1.set_crossings(l.link().data()[0..3].iter().cloned());
        b1.set_elements([]);
        b1.process_all();

        b2.set_crossings(l.link().data()[3..].iter().cloned());
        b2.set_elements([]);
        b2.process_all();

        b1.merge(b2);
        b1.deloop_all(false);
        b1.finalize();

        let c = b1.into_khi_complex();
        c.check_d_all();

        // TODO check homology
    }
}