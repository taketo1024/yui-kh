use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use num_traits::Zero;
use yui::bitseq::Bit;
use yui::{hashmap, Ring, RingOps};
use yui_link::{Crossing, Edge, Link};

use crate::ext::LinkExt;
use crate::kh::v2::cob::LcCobTrait;
use crate::kh::v2::tng::Tng;
use crate::kh::{KhAlgGen, KhChain, KhComplex, KhGen};

use super::cob::{Bottom, Dot, Cob, CobComp, LcCob};
use super::tng::TngComp;
use super::tng_complex::{TngComplex, TngKey};

pub struct TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    crossings: Vec<Crossing>,
    complex: TngComplex<R>,
    elements: Vec<Elem<R>>,
    pub auto_deloop: bool,
    pub auto_elim: bool
}

impl<R> From<TngComplex<R>> for TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(complex: TngComplex<R>) -> Self {
        Self { crossings: vec![], complex, elements: vec![], auto_deloop: true, auto_elim: true }
    }
}

impl<R> TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn build_kh_complex(l: &Link, h: &R, t: &R, reduced: bool) -> KhComplex<R> { 
        let deg_shift = KhComplex::deg_shift_for(l, reduced);
        let base_pt = if reduced { 
            l.first_edge()
        } else { 
            None
        };

        let mut b = Self::new(h, t, deg_shift, base_pt);

        if t.is_zero() && l.is_knot() {
            b.elements = Self::make_canon_cycles(l, base_pt);
        }
        
        b.set_crossings(l.data().clone());
        b.process_all();

        let canon_cycles = b.eval_elements();
        let complex = b.into_tng_complex().into_kh_complex(canon_cycles);

        complex
    }

    pub fn new(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let crossings = vec![];
        let complex = TngComplex::init(h, t, deg_shift, base_pt);
        let elements = vec![];

        let auto_deloop = true;
        let auto_elim   = true;

        Self { crossings, complex, elements, auto_deloop, auto_elim }
    }

    pub fn complex(&self) -> &TngComplex<R> { 
        &self.complex
    }

    pub fn set_crossings<I>(&mut self, crossings: I)
    where I: IntoIterator<Item = Crossing> {
        self.crossings = crossings.into_iter().collect_vec();
    }

    pub fn choose_next(&mut self) -> Option<Crossing> { 
        if self.crossings.is_empty() { 
            return None;
        }
        if self.crossings.len() == 1 { 
            let x = self.crossings.remove(0);
            return Some(x);
        }

        let mut candidate = 0;
        let mut score = 0;

        for (i, x) in self.crossings.iter().enumerate() { 
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

            let s = self.complex.iter_verts().map(|(_, v)| {
                v.tng().comps().map(|c| 
                    arcs.iter().filter(|a| 
                        c.path().is_connectable_bothends(a)
                    ).count()
                ).sum::<usize>()
            }).sum::<usize>();

            if s > score { 
                candidate = i;
                score = s;
            }
        }

        let x = self.crossings.remove(candidate);
        Some(x)
    }

    pub fn process_all(&mut self) { 
        while let Some(x) = self.choose_next() { 
            self.process(x)
        }
        self.finalize();
    }

    pub fn process(&mut self, x: Crossing) { 
        self.complex.append(&x);

        for e in self.elements.iter_mut() { 
            e.append(&x);
        }

        if self.auto_deloop { 
            while let Some((k, r)) = self.find_good_loop(false) { 
                self.deloop(&k, r);
            }
        }
    }

    pub fn find_loop(&self, allow_based: bool) -> Option<(TngKey, usize)> { 
        self.complex.iter_verts().find_map(|(k, v)| {
            v.tng().find_comp(|c| 
                c.is_circle() && (allow_based || !self.complex.contains_base_pt(c))
            ).map(|r| (*k, r))
        })
    }

    pub fn find_good_loop(&self, allow_based: bool) -> Option<(TngKey, usize)> { 
        let find_in = |k: &TngKey, loc_tng: &Tng| { 
            if let Some(r_loc) = loc_tng.find_comp(|c| { 
                c.is_circle() && (allow_based || !self.complex.contains_base_pt(c))
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
            v.out_edges().find_map(|l|
                self.complex.edge(k, l).gens().find_map(|cob| 
                    cob.comps().find_map(|c| 
                        if c.is_cap() || c.is_merge() { 
                            find_in(k, c.src())
                        } else if c.is_cup() || c.is_split() {
                            find_in(l, c.tgt())
                        } else { 
                            None
                        }
                    )
                )
            )
        )
    }

    pub fn deloop(&mut self, k: &TngKey, r: usize) {
        let c = self.complex.vertex(k).tng().comp(r);
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
        self.eliminate_elements(i, j);
        self.complex.eliminate(i, j);
    }

    fn eliminate_elements(&mut self, i: &TngKey, j: &TngKey) {
        let mut elements = std::mem::take(&mut self.elements);
        let a = self.complex.edge(i, j);
        for e in elements.iter_mut() { 
            let i_out = self.complex.keys_out_from(i).map(|k| 
                (k, self.complex.edge(i, k))
            );
            e.eliminate(i, j, a, i_out);
        }
        self.elements = elements;
    }

    pub fn finalize(&mut self) { 
        for b in [false, true] { 
            while let Some((k, r)) = self.find_loop(b) { 
                self.deloop(&k, r);
            }
        }    
        for e in self.elements.iter_mut() { 
            e.finalize();
        }
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }

    fn make_canon_cycles(l: &Link, base_pt: Option<Edge>) -> Vec<Elem<R>> { 
        assert_eq!(l.components().len(), 1);

        let p = base_pt.or(l.first_edge()).unwrap();
        let s = l.ori_pres_state().iter().enumerate().map(|(i, b)|
            (l.crossing_at(i).clone(), b)
        ).collect::<HashMap<_, _>>();

        let ori = if base_pt.is_some() { 
            vec![true]
        } else { 
            vec![true, false]
        };

        ori.into_iter().map(|o| { 
            let circles = l.colored_seifert_circles(p);
            let cob = Cob::new(
                circles.into_iter().map(|(circ, col)| { 
                    let t = TngComp::from(circ);
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
            Elem::new(cob, s.clone(), base_pt)
        }).collect()
    }

    pub fn eval_elements(&self) -> Vec<KhChain<R>> {
        let (h, t) = self.complex.ht();
        self.elements.iter().map(|z|
            z.eval(h, t, self.complex.deg_shift())
        ).collect()
    }
}

struct Elem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    init_cob: Cob,                     // precomposed at the final step.
    retr_cob: HashMap<TngKey, LcCob<R>>, // src must partially match init_cob. 
    state: HashMap<Crossing, Bit>,
    base_pt: Option<Edge>
}

impl<R> Elem<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(init_cob: Cob, state: HashMap<Crossing, Bit>, base_pt: Option<Edge>) -> Self { 
        let f = LcCob::from(Cob::empty());
        let retr_cob = hashmap! { TngKey::init() => f };
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

        let (k0, f0) = self.deloop_for(k, &f, c, KhAlgGen::X, Dot::None);
        self.retr_cob.insert(k0, f0);

        if !marked { 
            let (k1, f1) = self.deloop_for(k, &f, c, KhAlgGen::I, Dot::Y);
            self.retr_cob.insert(k1, f1);    
        }
    }

    fn deloop_for(&self, k: &TngKey, f: &LcCob<R>, c: &TngComp, label: KhAlgGen, dot: Dot) -> (TngKey, LcCob<R>) { 
        let mut k_new = *k;
        k_new.label.push(label);

        let f_new = f.clone().cap_off(Bottom::Tgt, c, dot);
        (k_new, f_new)
    }

    pub fn eliminate<'a, I>(&mut self, i: &TngKey, j: &TngKey, a: &LcCob<R>, i_out: I)
    where I: IntoIterator<Item = (&'a TngKey, &'a LcCob<R>)> {
        // mors into i can be simply dropped.
        self.retr_cob.remove(i);

        // mors into j must be redirected by -ca^{-1}
        let Some(f) = self.retr_cob.remove(j) else { return };
        let ainv = a.inv().unwrap();

        for (k, c) in i_out { 
            if k == j { continue }

            let caf = c * &ainv * &f;
            let d = if let Some(d) = self.retr_cob.remove(k) { 
                d - caf
            } else { 
                -caf
            };

            if !d.is_zero() { 
                self.retr_cob.insert(*k, d);
            }
        }
    }

    pub fn finalize(&mut self) { 
        let val = std::mem::take(&mut self.init_cob);
        let mut mors = std::mem::take(&mut self.retr_cob);

        mors = mors.into_iter().map(|(k, f)|
            (k, f.map_cob(|c| *c = &*c * &val))
        ).collect();
        mors.retain(|_, f| !f.is_zero());

        self.retr_cob = mors;
    }

    pub fn eval(&self, h: &R, t: &R, deg_shift: (isize, isize)) -> KhChain<R> {
        assert!(self.init_cob.is_empty());
        assert!(self.retr_cob.values().all(|f| f.is_closed()));

        KhChain::from_iter(self.retr_cob.iter().map(|(k, f)| {
            let x = KhGen::new(k.state, k.label, deg_shift);
            (x, f.eval(h, t))
        }))
    }
}

impl<R> Display for Elem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mors = self.retr_cob.iter().map(|(k, f)| { 
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
        let c = TngComplexBuilder::build_kh_complex(&l, &2, &0, false);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1() {
        let l = Link::from_pd_code([[0,0,1,1]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1_neg() {
        let l = Link::from_pd_code([[0,1,1,0]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::from_pd_code([[1,4,2,1],[2,4,3,3]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from_pd_code(pd_code);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 4);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_tangle() { 
        let mut c = TngComplexBuilder::new(&0, &0, (0, 0), None);
        c.set_crossings([
            Crossing::from_pd_code([4,2,5,1]),
            Crossing::from_pd_code([3,6,4,1])
        ]);

        c.process_all();
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::hopf_link();
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false);

        c.check_d_all();

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_8_19() {
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        let c = TngComplexBuilder::build_kh_complex(&l, &0, &0, false);

        c.check_d_all();

        let h = c.homology(false);

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
        let c = TngComplexBuilder::build_kh_complex(&l, &2, &0, false);
        let zs = c.canon_cycles();

        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);
        
        for z in zs { 
            assert!(z.gens().all(|x| x.h_deg() == 0));
            assert!(c.d(0, &z).is_zero());
        }
    }
}