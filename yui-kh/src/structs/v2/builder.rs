use std::collections::HashSet;

use log::info;
use yui::{Ring, RingOps};
use yui_link::{Link, Crossing, Edge};

use crate::ext::LinkExt;
use crate::{KhComplex, KhEnhState};

use super::tng::TngComp;
use super::cob::{Cob, CobComp, Dot};
use super::complex::TngComplex;
use super::elem::TngElem;

pub struct TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    crossings: Vec<Crossing>,
    complex: TngComplex<R>,
    canon_cycles: Vec<TngElem<R>>,
    base_pt: Option<Edge>,
    reduced: bool,
    pub auto_deloop: bool,
    pub auto_elim: bool
}

impl<R> TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        let base_pt = if reduced { 
            assert!(l.components().len() > 0);
            l.first_edge()
        } else { 
            None
        };

        let crossings = Self::sort_crossings(l, &base_pt);
        let deg_shift = KhComplex::<i64>::deg_shift_for(l, reduced);
        let complex = TngComplex::new(h, t, deg_shift);

        let auto_deloop = true;
        let auto_elim   = true;

        Self { crossings, complex, canon_cycles: vec![], base_pt, reduced, auto_deloop, auto_elim }
    }

    fn sort_crossings(l: &Link, base_pt: &Option<Edge>) -> Vec<Crossing> { 
        let mut remain = l.data().clone();
        let mut endpts = HashSet::new();
        let mut res = Vec::new();

        fn take_best(remain: &mut Vec<Crossing>, endpts: &mut HashSet<Edge>, base_pt: &Option<Edge>) -> Option<Crossing> { 
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
            res.push(x);
        }

        res
    }

    pub fn complex(&self) -> &TngComplex<R> { 
        &self.complex
    }

    pub fn into_complex(self) -> TngComplex<R> { 
        self.complex
    }

    pub fn canon_cycles(&self) -> &Vec<TngElem<R>> { 
        &self.canon_cycles
    }

    pub fn process(&mut self) {
        for i in 0 .. self.crossings.len() { 
            self.proceed_each(i);
        }
        self.finalize();
    }

    fn proceed_each(&mut self, i: usize) { 
        let x = &self.crossings[i];
        
        self.complex.append(x);

        for e in self.canon_cycles.iter_mut() { 
            e.append(x);
        }

        if self.auto_deloop { 
            let p = self.base_pt;
            while let Some((k, r, _)) = self.complex.find_loop(p) { 
                self.deloop(&k, r, false);
            }
        }
    }

    fn deloop(&mut self, k: &KhEnhState, r: usize, reduced: bool) {
        let c = self.complex.vertex(k).tng().comp(r);
        for e in self.canon_cycles.iter_mut() { 
            e.deloop(k, &c, reduced);
        }

        let keys = self.complex.deloop(k, r, reduced);

        if self.auto_elim { 
            for k in keys { 
                self.eliminate(&k)
            }
        }
    }

    fn eliminate(&mut self, k: &KhEnhState) {
        if let Some((i, j)) = self.complex.find_inv_edge(k) { 
            let i_out = self.complex.vertex(&i).out_edges();
            for e in self.canon_cycles.iter_mut() { 
                e.eliminate(&i, &j, i_out);
            }
            
            self.complex.eliminate(&i, &j);
        }
    }

    fn finalize(&mut self) { 
        if self.reduced && self.auto_deloop { 
            while let Some((k, r, _)) = self.complex.find_loop(None) { 
                self.deloop(&k, r, true);
            }
        }
        
        for e in self.canon_cycles.iter_mut() { 
            e.finalize();
        }
    }

    pub fn make_canon_cycles(&mut self) { 
        let l = Link::new(self.crossings.clone());
        
        assert_eq!(l.components().len(), 1);

        let s = l.ori_pres_state();
        let ori = if self.reduced { 
            vec![true]
        } else { 
            vec![true, false]
        };

        let cycles = ori.into_iter().map(|o| { 
            let circles = l.colored_seifert_circles(o);
            let f = Cob::new(
                circles.iter().map(|(circ, col)| { 
                    let t = TngComp::from(circ);
                    let mut cup = CobComp::cup(t);
                    let dot = if col.is_a() { Dot::X } else { Dot::Y };
                    cup.add_dot(dot);
                    cup
                }).collect()
            );
    
            let z = TngElem::init(s, f);
            info!("canon-cycle: {z}");
    
            z
        }).collect();
        
        self.canon_cycles = cycles;
    }
}

#[cfg(test)]
mod tests { 
    use yui_homology::{RModStr, ChainComplexCommon};
    use super::*;

    #[test]
    fn test_unknot_rm1() {
        let l = Link::from_pd_code([[0,0,1,1]]);
        let mut b = TngComplexBuilder::new(&l, &0, &0, false);
        b.process();

        let c = b.complex.eval(&0, &0);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1_neg() {
        let l = Link::from_pd_code([[0,1,1,0]]);
        let mut b = TngComplexBuilder::new(&l, &0, &0, false);
        b.process();

        let c = b.complex.eval(&0, &0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::from_pd_code([[1,4,2,1],[2,4,3,3]]);
        let mut b = TngComplexBuilder::new(&l, &0, &0, false);
        b.process();

        let c = b.complex.eval(&0, &0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from_pd_code(pd_code);
        let mut b = TngComplexBuilder::new(&l, &0, &0, false);
        b.process();

        let c = b.complex.eval(&0, &0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 4);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::hopf_link();
        let mut b = TngComplexBuilder::new(&l, &0, &0, false);
        b.process();

        let c = b.complex.eval(&0, &0);

        c.check_d_all();

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_8_19() {
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        let mut b = TngComplexBuilder::new(&l, &0, &0, false);
        b.process();

        let c = b.complex.eval(&0, &0);

        c.check_d_all();

        let h = c.homology(false);

        for i in [1,6,7,8] {
            assert_eq!(h[i].rank(), 0);
            assert_eq!(h[i].is_free(), true);
        }

        for i in [0,4,5] {
            assert_eq!(h[i].rank(), 2);
            assert_eq!(h[i].is_free(), true);
        }

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn canon_cycle_trefoil() { 
        let l = Link::trefoil();
        let mut b = TngComplexBuilder::new(&l, &0, &0, false);

        b.make_canon_cycles();
        b.process();

        assert_eq!(b.canon_cycles.len(), 2);
        
        for i in [0, 1] { 
            let z = &b.canon_cycles[i];
            let z = z.eval(&2, &0);
            println!("a[{i}] = {z}");
        }
    }
}