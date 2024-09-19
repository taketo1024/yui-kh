use std::collections::{HashMap, HashSet};
use std::fmt::Display;

use itertools::Itertools;
use num_traits::Zero;
use yui::bitseq::Bit;
use yui::{hashmap, Ring, RingOps};
use yui_homology::XChainComplex;
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
    complex: TngComplex<R>,
    crossings: Vec<Crossing>,
    elements: Vec<Elem<R>>,
    pub auto_deloop: bool,
    pub auto_elim: bool
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
        b.crossings = l.data().clone();

        if t.is_zero() && l.is_knot() {
            b.elements = Self::make_canon_cycles(l, base_pt);
        }
        
        b.process();

        let canon_cycles = b.eval_elements();
        let complex = b.into_raw_complex();

        KhComplex::new_impl(complex, canon_cycles, reduced, deg_shift)

    }

    pub fn new(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let complex = TngComplex::new(h, t, deg_shift, base_pt);
        let crossings = vec![];
        let elements = vec![];

        let auto_deloop = true;
        let auto_elim   = true;

        Self { complex, crossings, elements, auto_deloop, auto_elim }
    }

    pub fn add_crossing(&mut self, x: &Crossing) { 
        self.crossings.push(x.clone())
    }

    pub fn set_crossings<I>(&mut self, xs: I) 
    where I: IntoIterator<Item = Crossing> { 
        self.crossings = xs.into_iter().collect()
    }

    pub fn process(&mut self) {
        self.sort_crossings();

        for i in 0 .. self.crossings.len() { 
            self.proceed_each(i);
        }
        
        self.finalize();
    }

    fn sort_crossings(&mut self) { 
        let mut remain = std::mem::take(&mut self.crossings);
        let mut endpts = HashSet::new();
        let mut sorted = Vec::new();
        let base_pt = self.complex.base_pt();

        fn take_best(remain: &mut Vec<Crossing>, endpts: &mut HashSet<Edge>, base_pt: Option<Edge>) -> Option<Crossing> { 
            if remain.is_empty() { 
                return None 
            }

            let mut cand_i = 0;
            let mut cand_c = 0;

            for (i, x) in remain.iter().enumerate() { 
                if let Some(e) = base_pt { 
                    if x.edges().contains(&e) { 
                        continue
                    }
                }

                let c = x.edges().iter().filter(|e| 
                    endpts.contains(e)
                ).count();

                if c == 4 { 
                    let x = remain.remove(i);
                    return Some(x);
                } else if c > cand_c { 
                    cand_i = i;
                    cand_c = c;
                }
            }

            let x = remain.remove(cand_i);
            for e in x.edges() { 
                endpts.insert(*e);
            }
            
            Some(x)
        }

        while let Some(x) = take_best(&mut remain, &mut endpts, base_pt) { 
            sorted.push(x);
        }

        self.crossings = sorted;
    }

    fn proceed_each(&mut self, i: usize) { 
        let x = &self.crossings[i];
        
        self.complex.append(x);

        for e in self.elements.iter_mut() { 
            e.append(x);
        }

        if self.auto_deloop { 
            while let Some((k, r)) = self.complex.find_loop(false) { 
                self.deloop(&k, r);
            }
        }
    }

    fn deloop(&mut self, k: &TngKey, r: usize) {
        let c = self.complex.vertex(k).tng().comp(r);
        for e in self.elements.iter_mut() { 
            e.deloop(k, c);
        }

        let keys = self.complex.deloop(k, r);

        if self.auto_elim { 
            for k in keys { 
                self.eliminate(&k)
            }
        }
    }

    fn eliminate(&mut self, k: &TngKey) {
        if let Some((i, j)) = self.complex.find_inv_edge(k) { 
            let i_out = self.complex.vertex(&i).out_edges();
            for e in self.elements.iter_mut() { 
                e.eliminate(&i, &j, i_out);
            }
            
            self.complex.eliminate(&i, &j);
        }
    }

    fn finalize(&mut self) { 
        if self.auto_deloop { 
            while let Some((k, r)) = self.complex.find_loop(true) { 
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

    pub fn into_raw_complex(self) -> XChainComplex<KhGen, R> { 
        self.complex.into_complex()
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
        let tng = Tng::from_resolved(&a, self.base_pt);
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

        let tng = Tng::from_resolved(x, self.base_pt);
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

    pub fn eliminate(&mut self, i: &TngKey, j: &TngKey, i_out: &HashMap<TngKey, LcCob<R>>) {
        // mors into i can be simply dropped.
        self.retr_cob.remove(i);

        // mors into j must be redirected by -ca^{-1}
        let Some(f) = self.retr_cob.remove(j) else { return };

        let a = &i_out[&j];
        let ainv = a.inv().unwrap();

        for (k, c) in i_out.iter() { 
            if k == j { continue }

            let caf = c * &ainv * &f;
            let d = if let Some(d) = self.retr_cob.get(k) { 
                d - caf
            } else { 
                -caf
            };

            if !d.is_zero() { 
                self.retr_cob.insert(*k, d);
            } else { 
                self.retr_cob.remove(k);
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