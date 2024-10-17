use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use log::{debug, info};
use num_traits::Zero;
use yui::bitseq::Bit;
use yui::{hashmap, Ring, RingOps};
use yui_link::{Crossing, Edge, Link};

use crate::ext::LinkExt;
use crate::kh::{KhAlgGen, KhChain, KhComplex, KhGen};

use super::cob::{Bottom, Dot, Cob, CobComp, LcCobTrait, LcCob};
use super::tng::{Tng, TngComp};
use super::tng_complex::{TngComplex, TngKey};

pub struct TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    crossings: Vec<Crossing>,
    complex: TngComplex<R>,
    elements: Vec<BuildElem<R>>,
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
        assert!(!reduced || l.components().len() > 0);

        let deg_shift = KhComplex::deg_shift_for(l, reduced);
        let base_pt = reduced.then(|| l.first_edge().unwrap());

        let mut b = Self::new(h, t, deg_shift, base_pt);
        b.set_crossings(l.data().clone());

        if t.is_zero() && l.is_knot() {
            let canon = Self::make_canon_cycles(l, base_pt, reduced);
            b.set_elements(canon);
        }

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

    pub(crate) fn set_elements<I>(&mut self, elements: I)
    where I: IntoIterator<Item = BuildElem<R>> { 
        self.elements = elements.into_iter().collect_vec();
    }

    pub fn choose_next(&mut self) -> Option<Crossing> { 
        let Some((i, _)) = self.crossings.iter().enumerate().max_by_key(|(_, x)|
            count_loops(&self.complex, x)
        ) else { 
            return None
        };

        let x = self.crossings.remove(i);
        Some(x)
    }

    pub fn process_all(&mut self) { 
        while let Some(x) = self.choose_next() { 
            self.process(&x)
        }
        self.finalize();
    }

    pub fn process(&mut self, x: &Crossing) { 
        info!("(n: {}, v: {}) append: {x}", self.complex.dim(), self.complex.nverts());

        self.complex.append(x);

        for e in self.elements.iter_mut() { 
            e.append(x);
        }

        if self.auto_deloop { 
            while let Some((k, r)) = self.find_loop(false) { 
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

    pub fn deloop(&mut self, k: &TngKey, r: usize) {
        let c = self.complex.vertex(k).tng().comp(r);

        info!("(n: {}, v: {}) deloop: {c} in {}", self.complex.dim(), self.complex.nverts(), self.complex.vertex(k));

        for (idx, e) in self.elements.iter_mut().enumerate() { 
            debug!("  e[{idx}] {e}");

            e.deloop(k, c);

            debug!("    -> {e}");
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
        info!("(n: {}, v: {}) eliminate {}: {} -> {}", self.complex.dim(), self.complex.nverts(), self.complex.edge(i, j), self.complex.vertex(i), self.complex.vertex(j));
        
        self.eliminate_elements(i, j);
        self.complex.eliminate(i, j);
    }

    fn eliminate_elements(&mut self, i: &TngKey, j: &TngKey) {
        let a = self.complex.edge(i, j);
        for (idx, e) in self.elements.iter_mut().enumerate() { 
            debug!("  e[{idx}] {e}");

            let i_out = self.complex.keys_out_from(i).map(|k| 
                (k, self.complex.edge(i, k))
            );
            e.eliminate(i, j, a, i_out);

            debug!("    -> {e}");
        }
    }

    pub fn finalize(&mut self) { 
        for b in [false, true] { 
            while let Some((k, r)) = self.find_loop(b) { 
                self.deloop(&k, r);
            }
        }    
    }

    pub(crate) fn take_elements(&mut self) -> Vec<BuildElem<R>> {
        std::mem::take(&mut self.elements)
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }

    pub(crate) fn make_canon_cycles(l: &Link, base_pt: Option<Edge>, reduced: bool) -> Vec<BuildElem<R>> { 
        assert!(l.is_knot());
        assert!(!reduced || base_pt.is_some());

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
}

pub(crate) fn count_loops<R>(complex: &TngComplex<R>, x: &Crossing) -> usize
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let arcs = if x.is_resolved() { 
        let a = x.arcs();
        vec![a.0, a.1]
    } else { 
        let a0 = x.resolved(Bit::Bit0).arcs();
        let a1 = x.resolved(Bit::Bit1).arcs();
        vec![a0.0, a0.1, a1.0, a1.1]
    }.into_iter().filter(|a|
        complex.base_pt().map(|e| !a.contains(e)).unwrap_or(true)
    ).collect_vec();

    let count = complex.iter_verts().map(|(_, v)| {
        v.tng().comps().map(|c| 
            arcs.iter().filter(|a| 
                c.path().is_connectable_bothends(a)
            ).count()
        ).sum::<usize>()
    }).sum::<usize>();

    count
}
pub(crate) struct BuildElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    init_cob: Cob,                       // initial cob, precomposed at the final step.
    retr_cob: HashMap<TngKey, LcCob<R>>, // building cob, src must always match init_cob. 
    state: HashMap<Crossing, Bit>,
    base_pt: Option<Edge>
}

impl<R> BuildElem<R> 
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

    pub fn eval(&self, h: &R, t: &R, deg_shift: (isize, isize)) -> KhChain<R> {
        let init = LcCob::from(self.init_cob.clone());

        debug!("  eval {init} : {}", self);

        assert!(self.retr_cob.values().all(|retr| init.is_stackable(retr)));

        let eval = self.retr_cob.iter().filter_map(|(k, retr)| {
            let f = retr * &init;
            (!f.is_zero()).then(|| (*k, f))
        }).map(|(k, f)| { 
            let x = KhGen::new(k.state, k.label, deg_shift);
            (x, f.eval(h, t))
        }).collect::<KhChain<R>>();

        debug!("    -> {eval}");

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
    use yui_homology::{ChainComplexTrait, SummandTrait};

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

        let h = c.inner().homology();

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