use std::collections::{BTreeSet, HashMap};
use std::fmt::Display;
use std::ops::RangeInclusive;

use itertools::Itertools;
use log::info;
use num_traits::Zero;
use yui::bitseq::Bit;
use yui::{hashmap, Ring, RingOps};
use yui_link::{Crossing, Edge, Link};

use crate::ext::LinkExt;
use crate::kh::{KhAlgGen, KhChain, KhComplex};

use super::cob::{Bottom, Dot, Cob, CobComp, LcCobTrait, LcCob};
use super::tng::{Tng, TngComp};
use super::tng_complex::{TngComplex, TngKey};

pub struct TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    crossings: Vec<Crossing>,
    complex: TngComplex<R>,
    elements: Vec<BuildElem<R>>,
    max_dim: usize,
    h_range: Option<RangeInclusive<isize>>,
    pub auto_deloop: bool,
    pub auto_elim: bool
}

impl<R> From<TngComplex<R>> for TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(complex: TngComplex<R>) -> Self {
        let max_dim = complex.dim();
        Self { 
            crossings: vec![], 
            complex, 
            elements: vec![], 
            max_dim, 
            h_range: None, 
            auto_deloop: true, 
            auto_elim: true 
        }
    }
}

impl<R> TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn build_kh_complex(l: &Link, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> KhComplex<R> { 
        let base_pt = if reduced { l.first_edge() } else { None };
        let mut b = Self::new(l, h, t, base_pt);
        b.h_range = h_range;
        b.process_all();
        b.finalize();
        b.into_kh_complex()
    }

    pub fn new(l: &Link, h: &R, t: &R, base_pt: Option<Edge>) -> Self { 
        let reduced = base_pt.is_some();
        let deg_shift = KhComplex::deg_shift_for(l, reduced);

        let mut b = Self::init(h, t, deg_shift, base_pt);
        b.set_crossings(l.data().clone());
        b.max_dim = l.crossing_num();

        if t.is_zero() && l.is_knot() {
            let canon = Self::make_canon_cycles(l, base_pt);
            b.set_elements(canon);
        }

        b
    }

    pub fn init(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let complex = TngComplex::init(h, t, deg_shift, base_pt);
        Self::from(complex)
    }

    pub fn complex(&self) -> &TngComplex<R> { 
        &self.complex
    }

    pub(crate) fn complex_mut(&mut self) -> &mut TngComplex<R> { 
        &mut self.complex
    }

    pub fn crossings(&self) -> impl Iterator<Item = &Crossing> { 
        self.crossings.iter()
    }

    pub fn set_crossings<I>(&mut self, crossings: I)
    where I: IntoIterator<Item = Crossing> {
        self.crossings = crossings.into_iter().collect_vec();
    }

    pub(crate) fn take_crossings(&mut self) -> Vec<Crossing> { 
        std::mem::take(&mut self.crossings)
    }

    pub fn elements(&self) -> impl Iterator<Item = &BuildElem<R>> { 
        self.elements.iter()
    }

    pub fn set_elements<I>(&mut self, elements: I)
    where I: IntoIterator<Item = BuildElem<R>> { 
        self.elements = elements.into_iter().collect_vec();
    }

    pub(crate) fn take_elements(&mut self) -> Vec<BuildElem<R>> {
        std::mem::take(&mut self.elements)
    }

    pub fn h_range(&self) -> Option<RangeInclusive<isize>> { 
        self.h_range.clone()
    }

    pub fn set_h_range(&mut self, h_range: RangeInclusive<isize>) { 
        self.h_range = Some(h_range)
    }

    pub fn choose_next(&mut self) -> Option<Crossing> { 
        let Some((i, _)) = self.crossings.iter().enumerate().max_by_key(|(_, x)|
            self.count_connections(x)
        ) else { 
            return None
        };

        let x = self.crossings.remove(i);
        Some(x)
    }

    fn count_connections(&self, x: &Crossing) -> usize { 
        let arcs = if x.is_resolved() { 
            let a = x.arcs();
            vec![a.0, a.1]
        } else { 
            let a0 = x.resolved(Bit::Bit0).arcs();
            let a1 = x.resolved(Bit::Bit1).arcs();
            vec![a0.0, a0.1, a1.0, a1.1]
        }.into_iter().filter(|a|
            self.complex.base_pt().map(|e| !a.contains(e)).unwrap_or(true)
        ).collect_vec();

        let count = self.complex.iter_verts().map(|(_, v)| {
            v.tng().comps().map(|c| 
                arcs.iter().map(|a| 
                    match c.path() {
                        p if p.is_connectable_bothends(a) => 2,
                        p if p.is_connectable(a)          => 1,
                        _                                 => 0
                    }
                ).sum::<usize>()
            ).sum::<usize>()
        }).sum::<usize>();

        count
    }

    pub fn process_all(&mut self) { 
        while let Some(x) = self.choose_next() { 
            self.process(&x)
        }
    }

    pub fn process(&mut self, x: &Crossing) { 
        info!("({}) append: {x}", self.stat());

        if let Some(i) = self.crossings.iter().find_position(|&e| e == x) { 
            self.crossings.remove(i.0);
        }
        self.complex.append(x);

        for e in self.elements.iter_mut() { 
            e.append(x);
        }

        if self.h_range.is_some() { 
            self.drop_vertices();
        }

        if self.auto_deloop { 
            self.deloop_all(false);
        }
    }

    pub fn drop_vertices(&mut self) { 
        let Some(h_range) = &self.h_range else { return };

        let i0 = self.complex.deg_shift().0;
        let remain = self.max_dim - self.complex.dim();

        let drop = self.complex.keys().filter(|k| 
            i0 + ((k.weight() + remain) as isize) < *h_range.start() || 
            i0 + (k.weight() as isize) > *h_range.end()
        ).cloned().collect_vec();

        for k in drop { 
            info!("({}) drop {}", self.stat(), self.complex.vertex(&k));
            self.complex.remove_vertex(&k);
        }

        // TODO drop elements
    }

    pub fn deloop_all(&mut self, allow_based: bool) { 
        let mut keys = self.complex.keys().filter(|k| 
            self.complex.vertex(k).tng().contains_circle()
        ).cloned().collect::<BTreeSet<_>>();

        for special in [true, false] { 
            while let Some((k, r)) = self.find_loop(allow_based, special, keys.iter()) { 
                keys.remove(&k);
    
                let updated = self.deloop(&k, r);
                
                keys.extend(updated);
                keys.retain(|k| 
                    self.complex.contains_key(k) && 
                    self.complex.vertex(k).tng().contains_circle()
                );
            }
        }
    }

    pub fn find_loop<'a, I>(&self, allow_based: bool, special: bool, keys: I) -> Option<(TngKey, usize)>
    where I: IntoIterator<Item = &'a TngKey> { 
        let ok = |c: &TngComp| { 
            c.is_circle() && (allow_based || !self.complex.contains_base_pt(c))
        };

        let plain_merge = |k: &TngKey, c: &TngComp| { 
            self.complex.keys_out_from(k).find(|l|
                self.complex.edge(k, l).iter().find(|(cob, a)|
                    a.is_unit() && cob.comps().find(|cob|
                        cob.src().contains(c) && cob.is_plain() && cob.is_merge()
                    ).is_some()
                ).is_some()
            ).is_some()
        };

        let plain_split = |k: &TngKey, c: &TngComp| {
            self.complex.keys_into(k).find(|j|
                self.complex.edge(j, k).iter().find(|(cob, a)|
                    a.is_unit() && cob.comps().find(|cob|
                        cob.tgt().contains(c) && cob.is_plain() && cob.is_split()
                    ).is_some()
                ).is_some()
            ).is_some()
        };

        let special_ok = |k: &TngKey, c: &TngComp| {
            ok(c) && (plain_merge(k, c) || plain_split(k, c))
        };

        keys.into_iter().find_map(|k|
            self.complex.vertex(k).tng().find_comp(|c|
                if special { 
                    special_ok(k, c)
                } else { 
                    ok(c)
                }
            ).map(|r| (*k, r))
        )
    }

    pub fn deloop(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex.vertex(k).tng().comp(r);

        info!("({}) deloop: {c} in {}", self.stat(), self.complex.vertex(k));

        for e in self.elements.iter_mut() { 
            e.deloop(k, c);
        }

        let keys = self.complex.deloop(k, r);

        if self.auto_elim { 
            for k in keys.iter() { 
                if let Some((i, j)) = self.find_inv_edge(k) { 
                    self.eliminate(&i, &j)
                }
            }
        }

        keys
    }

    pub fn find_inv_edge(&self, k: &TngKey) -> Option<(TngKey, TngKey)> { 
        let v = self.complex.vertex(k);
        let edges = Iterator::chain(
            v.in_edges().map(|i| (i, k)),
            v.out_edges().map(|i| (k, i))
        );

        edges.filter_map(|(i, j)| {
            let f = self.complex.edge(i, j);
            f.is_invertible().then(|| {
                let s = self.complex.elim_weight(i, j);
                (s, i, j)
            })
        })
        .min_by_key(|(s, _, _)| *s)
        .map(|(_, i, j)| (*i, *j))
    }

    pub fn eliminate(&mut self, i: &TngKey, j: &TngKey) {
        info!("({}) eliminate {}: {} -> {}", self.stat(), self.complex.edge(i, j), self.complex.vertex(i), self.complex.vertex(j));
        
        self.eliminate_elements(i, j);
        self.complex.eliminate(i, j);
    }

    fn eliminate_elements(&mut self, i: &TngKey, j: &TngKey) {
        let a = self.complex.edge(i, j);
        for e in self.elements.iter_mut() { 
            let i_out = self.complex.keys_out_from(i).map(|k| 
                (k, self.complex.edge(i, k))
            );
            e.eliminate(i, j, a, i_out);
        }
    }

    pub fn finalize(&mut self) { 
        info!("finalize");

        self.deloop_all(false);
        self.deloop_all(true);

        assert!(self.complex.is_completely_delooped());
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }

    pub fn into_kh_complex(self) -> KhComplex<R> { 
        let h_range = self.h_range.clone();
        let canon_cycles = self.eval_elements();
        let ckh = self.into_tng_complex().into_kh_complex(canon_cycles);

        if let Some(h_range) = h_range { 
            ckh.truncated(h_range)
        } else { 
            ckh
        }
    }

    pub(crate) fn make_canon_cycles(l: &Link, base_pt: Option<Edge>) -> Vec<BuildElem<R>> { 
        assert!(l.is_knot());

        let reduced = base_pt.is_some();
        let start_p = base_pt.or(l.first_edge()).unwrap();
        let circles = l.colored_seifert_circles(start_p);

        let state_map = l.ori_pres_state().iter().enumerate().map(|(i, b)|
            (l.crossing_at(i).clone(), b)
        ).collect::<HashMap<_, _>>();

        let ori = if reduced { 
            vec![true]
        } else { 
            vec![true, false]
        };

        let cycles = ori.into_iter().map(|o| { 
            let cob = Cob::new(
                circles.iter().map(|(circ, col)| { 
                    let t = TngComp::from(circ.clone());
                    let mut cup = CobComp::cup(t);
                    let dot = if col.is_a() == o { 
                        Dot::X 
                    } else { 
                        Dot::Y 
                    };
                    cup.add_dot(dot);
                    cup
                })
            );
            BuildElem::new(cob, state_map.clone(), base_pt)
        }).collect();

        cycles
    }

    pub fn eval_elements(&self) -> Vec<KhChain<R>> {
        let (h, t) = self.complex.ht();
        self.elements.iter().map(|z|
            z.eval(h, t, self.complex.deg_shift())
        ).collect()
    }

    pub(crate) fn stat(&self) -> String { 
        self.complex.stat()
    }
}

#[derive(Clone)]
pub struct BuildElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    init_cob: Cob,                       // initial cob, precomposed at the final step.
    retr_cob: HashMap<TngKey, LcCob<R>>, // building cob, src must always match init_cob. 
    state: HashMap<Crossing, Bit>,
    base_pt: Option<Edge>
}

impl<R> BuildElem<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(init_cob: Cob, state: HashMap<Crossing, Bit>, base_pt: Option<Edge>) -> Self { 
        let k0 = TngKey::init();
        let f0 = LcCob::from(Cob::empty());
        let retr_cob = hashmap! { k0 => f0 };
        Self{ init_cob, retr_cob, state, base_pt }
    }

    pub fn append(&mut self, x: &Crossing) { 
        if !x.is_resolved() { 
            self.append_x(x)
        } else { 
            self.append_a(x)
        }
    }

    fn append_x(&mut self, x: &Crossing) {
        assert!(!x.is_resolved());

        let r = self.state[x];
        let a = x.resolved(r);
        let tng = Tng::from_resolved(&a);
        let id = Cob::id(&tng);

        let mors = std::mem::take(&mut self.retr_cob);
        self.retr_cob = mors.into_iter().map(|(mut k, f)| {
            k.state.push(r);
            let f = f.connect(&id);
            (k, f)
        }).collect();
    }

    fn append_a(&mut self, x: &Crossing) {
        assert!(x.is_resolved());

        let tng = Tng::from_resolved(x);
        let id = Cob::id(&tng);

        let mors = std::mem::take(&mut self.retr_cob);
        self.retr_cob = mors.into_iter().map(|(k, f)| {
            let f = f.connect(&id);
            (k, f)
        }).collect();
    }

    pub fn deloop(&mut self, k: &TngKey, c: &TngComp) {
        let Some(f) = self.retr_cob.remove(k) else { return };
        let marked = self.base_pt.map(|e| c.contains(e)).unwrap_or(false);

        let k0 = k + KhAlgGen::X;
        let f0 = f.clone().cap_off(Bottom::Tgt, c, Dot::None);
        self.retr_cob.insert(k0, f0);

        if !marked { 
            let k1 = k + KhAlgGen::I;
            let f1 = f.cap_off(Bottom::Tgt, c, Dot::Y);
            self.retr_cob.insert(k1, f1);    
        }
    }

    //  Gaussian Elimination
    //
    //       a
    //  v0 - - -> v1         .             .
    //     \   / b
    //       /         ==>  
    //     /   \ c              d - ca⁻¹b
    //  w0 -----> w1         w0 ---------> w1
    //       d                
    
    pub fn eliminate<'a, I>(&mut self, i: &TngKey, j: &TngKey, a: &LcCob<R>, i_out: I)
    where I: IntoIterator<Item = (&'a TngKey, &'a LcCob<R>)> {
        // mors into i can be simply dropped.
        self.retr_cob.remove(i);

        // mors into j must be redirected by -ca^{-1}
        let Some(b) = self.retr_cob.remove(j) else { return };
        let ainv = a.inv().unwrap();

        for (k, c) in i_out { 
            if k == j { continue }

            let cab = c * &ainv * &b;
            let s = if let Some(d) = self.retr_cob.remove(k) { 
                d - cab
            } else { 
                -cab
            };

            if !s.is_zero() { 
                self.retr_cob.insert(*k, s);
            }
        }
    }

    pub fn is_evalable(&self) -> bool { 
        let init = LcCob::from(self.init_cob.clone());
        self.retr_cob.values().all(|c| init.is_stackable(c)) && 
        self.retr_cob.values().all(|c| c.iter().all(|(c, _)| c.tgt().is_empty()))
    }

    pub fn eval(&self, h: &R, t: &R, deg_shift: (isize, isize)) -> KhChain<R> {
        assert!(self.is_evalable());

        let init = LcCob::from(self.init_cob.clone());
        let eval = self.retr_cob.iter().map(|(k, retr)| {
            let x = k.as_gen(deg_shift);
            let f = retr * &init;
            let r = f.eval(h, t);
            (x, r)
        }).collect::<KhChain<R>>();

        eval
    }

    pub fn modify<F>(&mut self, f: F)
    where F: Fn(TngKey, LcCob<R>) -> (TngKey, LcCob<R>) { 
        let retr_cob = std::mem::take(&mut self.retr_cob);
        self.retr_cob = retr_cob.into_iter().map(|(k, cob)|
            f(k, cob)
        ).collect();
    }
}

impl<R> Display for BuildElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mors = self.retr_cob.iter().sorted_by_key(|(&k, _)| k).map(|(k, f)| { 
            format!("{}: {}", k, f)
        }).join(", ");
        write!(f, "[{}]", mors)
    }
}

#[cfg(test)]
mod tests { 
    use num_traits::Zero;
    use yui_homology::{ChainComplexCommon, RModStr};

    use super::*;

    #[test]
    fn test_unknot() {
        let l = Link::unknot();
        let c = TngComplexBuilder::build_kh_complex(&l, &2, &0, false, None);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1() {
        let l = Link::from_pd_code([[0,0,1,1]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false, None);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1_neg() {
        let l = Link::from_pd_code([[0,1,1,0]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false, None);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::from_pd_code([[1,4,2,1],[2,4,3,3]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false, None);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from_pd_code(pd_code);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false, None);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 4);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_tangle() { 
        let mut c = TngComplexBuilder::init(&0, &0, (0, 0), None);
        c.set_crossings([
            Crossing::from_pd_code([4,2,5,1]),
            Crossing::from_pd_code([3,6,4,1])
        ]);

        c.process_all();
        
        assert!(!c.complex.is_completely_delooped());
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::hopf_link();
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false, None);

        c.check_d_all();

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_8_19() {
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false, None);

        c.check_d_all();

        let h = c.inner().homology(false);

        for i in [1,6,7,8] {
            assert_eq!(h[i].rank(), 0);
            assert!(h[i].is_free());
        }

        for i in [0,4,5] {
            assert_eq!(h[i].rank(), 2);
            assert!(h[i].is_free());
        }

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn canon_cycle_trefoil() { 
        let l = Link::trefoil();
        let c = TngComplexBuilder::build_kh_complex(&l, &2, &0, false, None);
        let zs = c.canon_cycles();

        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);
        
        for z in zs { 
            assert!(z.gens().all(|x| x.h_deg() == 0));
            assert!(c.d(0, &z).is_zero());
        }
    }

    #[test]
    fn partial_cpx() { 
        let l = Link::load("4_1").unwrap();
        let range = -1..=1;

        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false, Some(range.clone()));
        assert_eq!(c.h_range(), range);
        c.check_d_all();

        let h = c.homology();
        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());
    }
}