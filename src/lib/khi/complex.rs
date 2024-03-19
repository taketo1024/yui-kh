use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use yui::bitseq::{Bit, BitSeq};
use yui::lc::Lc;
use yui::{EucRing, EucRingOps, Ring, RingOps};
use yui_homology::{Grid, XChainComplex, XHomology, XModStr};
use yui_link::{Edge, State};

use crate::v1::cube::KhCube;
use crate::{KhComplex, KhGen, KhLabel};
use super::{InvLink, KhIGen};

pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    cube: KhCube<R>,
    state_map: HashMap<State, State>,
    label_map: HashMap<State, HashMap<usize, usize>>,
    deg_shift: (isize, isize)
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R, reduce_e: Option<Edge>) -> Self { 
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");

        let n = l.link().crossing_num();
        let deg_shift = KhComplex::deg_shift_for(l.link(), false);

        let t = R::zero();
        let cube = KhCube::new(l.link(), h, &t, reduce_e);

        let state_map = State::generate(n).map(|s| { 
            let xs = s.iter().enumerate().filter_map(|(i, b)| 
                b.is_one().then_some(i)
            ).map(|i| 
                l.inv_x(i)
            ).collect::<HashSet<_>>();

            let t = State::from_iter((0..n).map(|i| 
                if xs.contains(&i) { Bit::Bit1 } else { Bit::Bit0 }
            ));
            (s, t)
        }).collect::<HashMap<_, _>>();

        let label_map = State::generate(n).map(|s| { 
            let t = state_map.get(&s).cloned().unwrap();
            let v = cube.vertex(&s);
            let w = cube.vertex(&t);

            debug_assert_eq!(v.circles().len(), w.circles().len());

            let r = v.circles().len();
            let map = (0..r).map(|i| {
                let c0 = &v.circles()[i];
                let e = l.inv_e(c0.min_edge());
                let Some((j, c1)) = w.circles().iter().find_position(|c| c.contains(e)) else { 
                    panic!()
                };

                debug_assert_eq!(c0.edges().len(), c1.edges().len());

                (i, j)
            }).collect::<HashMap<_, _>>();

            (s, map)
        }).collect::<HashMap<_, _>>();

        Self { cube, state_map, label_map, deg_shift }
    }

    fn t_state(&self, s: State) -> State { 
        self.state_map.get(&s).cloned().unwrap()
    }

    fn t_label(&self, s: State, l: KhLabel) -> KhLabel { 
        debug_assert_eq!(l.len(), self.cube.vertex(&s).circles().len());

        let map = self.label_map.get(&s).unwrap();
        let mut seq = BitSeq::zeros(l.len());

        for (i, e) in  l.iter().enumerate() { 
            if e.is_1() { 
                let j = map.get(&i).cloned().unwrap();
                seq.set_1(j);
            }
        }

        KhLabel(seq)
    }

    fn t(&self, x: &KhGen) -> KhGen {
        let s = x.state;
        let t = self.t_state(s);
        let l = self.t_label(s, x.label);
        KhGen::new(t, l)
    }

    // f = 1 + τ
    fn f(&self, x: &KhGen) -> Lc<KhIGen, R> { 
        let qx = KhIGen::Q(x.clone());
        let qtx = KhIGen::Q(self.t(x));
        Lc::from_iter([
            (qx,  R::one()), 
            (qtx, R::one())
        ])
    }

    pub fn generators(&self, i: isize) -> Vec<KhIGen> { 
        let i0 = self.deg_shift.0;
        let i = i - i0;

        let b_gens = self.cube.generators(i).iter().map(|&&x| 
            KhIGen::B(x)
        ).collect_vec();

        let q_gens = if i > 0 { 
            self.cube.generators(i - 1).iter().map(|&&x|
                KhIGen::Q(x)
            ).collect()
        } else { 
            vec![]
        };

        let mut gens = b_gens;
        gens.extend(q_gens);
        gens
    }

    pub fn summand(&self, i: isize) -> XModStr<KhIGen, R> { 
        XModStr::free(self.generators(i))
    }

    pub fn differentiate(&self, z: &Lc<KhIGen, R>) -> Lc<KhIGen, R> { 
        z.apply(|x| self.d(x))
    }

    fn d(&self, x: &KhIGen) -> Lc<KhIGen, R> { 
        match x {
            KhIGen::B(x) => {
                let dx = self.cube.d(x).map_gens(|&y| KhIGen::B(y));
                let fx = self.f(x);
                dx + fx
            },
            KhIGen::Q(x) => {
                self.cube.d(x).map_gens(|&y| KhIGen::Q(y))
            },
        }
    }

    pub fn into_complex(self) -> XChainComplex<KhIGen, R> {
        let i0 = self.deg_shift.0;
        let range = self.cube.h_range();
        let range = (range.start() + i0) ..= (range.end() + i0 + 1);

        XChainComplex::new(
            Grid::generate(range, |i| self.summand(i)),
            1, 
            move |_, z| self.differentiate(z)
        )
    }   
}

impl<R> KhIComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(self, with_trans: bool) -> XHomology<KhIGen, R> {
        self.into_complex().reduced().homology(with_trans)
    }
}


#[cfg(test)]
mod tests {
    #![allow(unused)]

    use yui::poly::Poly;
    use yui::FF;
    use num_traits::{Zero, One};
    use yui_homology::{ChainComplexCommon, ChainComplexTrait, DisplaySeq};
    use yui_link::Link;
    use crate::KhHomology;

    use super::*;
 
    #[test]
    fn tau_state() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None);

        assert_eq!(
            c.t_state(State::from([0,0,0])), 
            State::from([0,0,0])
        );
        
        assert_eq!(
            c.t_state(State::from([1,0,0])), 
            State::from([1,0,0])
        );

        assert_eq!(
            c.t_state(State::from([0,1,0])), 
            State::from([0,0,1])
        );

        assert_eq!(
            c.t_state(State::from([0,0,1])), 
            State::from([0,1,0])
        );

        assert_eq!(
            c.t_state(State::from([1,1,0])), 
            State::from([1,0,1])
        );

        assert_eq!(
            c.t_state(State::from([1,0,1])), 
            State::from([1,1,0])
        );

        assert_eq!(
            c.t_state(State::from([0,1,1])), 
            State::from([0,1,1])
        );

        assert_eq!(
            c.t_state(State::from([1,1,1])), 
            State::from([1,1,1])
        );
    }

    #[test]
    fn tau_label() { 
        use crate::KhAlgGen::*;

        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None);

        assert_eq!(
            c.t_label(State::from([0,0,0]), KhLabel::from([I, X])), 
            KhLabel::from([I, X])
        );

        assert_eq!(
            c.t_label(State::from([0,1,0]), KhLabel::from([X])), 
            KhLabel::from([X])
        );

        assert_eq!(
            c.t_label(State::from([1,1,0]), KhLabel::from([I, X])), 
            KhLabel::from([X, I])
        );

        assert_eq!(
            c.t_label(State::from([1,1,1]), KhLabel::from([I, I, X])), 
            KhLabel::from([X, I, I])
        );
    }

    #[test]
    fn generators() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None);

        assert_eq!(c.generators(0).len(), 4);
        assert_eq!(c.generators(1).len(), 10);
        assert_eq!(c.generators(2).len(), 18);
        assert_eq!(c.generators(3).len(), 20);
        assert_eq!(c.generators(4).len(), 8);
    }

    #[test]
    fn d_lower_sym() { 
        use crate::KhAlgGen::*;

        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None);

        let x = KhIGen::B(
            KhGen::new(
                State::from([0,0,0]),
                KhLabel::from([X, I])
            )
        );
        let dx = c.d(&x);

        assert_eq!(dx, Lc::from_iter([
            (KhIGen::B(
                KhGen::new(
                    State::from([1,0,0]),
                    KhLabel::from([X])
                )
            ), R::one()),

            (KhIGen::B(
                KhGen::new(
                    State::from([0,1,0]),
                    KhLabel::from([X])
                )
            ), R::one()),
            
            (KhIGen::B(
                KhGen::new(
                    State::from([0,0,1]),
                    KhLabel::from([X])
                )
            ), R::one())
        ]));
    }

    #[test]
    fn d_lower_asym() { 
        use crate::KhAlgGen::*;

        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None);

        let x = KhIGen::B(
            KhGen::new(
                State::from([0,1,0]),
                KhLabel::from([X])
            )
        );
        let dx = c.d(&x);

        assert_eq!(dx, Lc::from_iter([
            (KhIGen::B(
                KhGen::new(
                    State::from([1,1,0]),
                    KhLabel::from([X,X])
                )
            ), R::one()),

            (KhIGen::B(
                KhGen::new(
                    State::from([0,1,1]),
                    KhLabel::from([X,X])
                )
            ), R::one()),
            
            (KhIGen::Q(
                KhGen::new(
                    State::from([0,1,0]),
                    KhLabel::from([X])
                )
            ), R::one()),

                        
            (KhIGen::Q(
                KhGen::new(
                    State::from([0,0,1]),
                    KhLabel::from([X])
                )
            ), R::one())
        ]));
    
    }

    #[test]
    fn d_upper() { 
        use crate::KhAlgGen::*;

        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None);

        let x = KhIGen::Q(
            KhGen::new(
                State::from([0,1,0]),
                KhLabel::from([X])
            )
        );
        let dx = c.d(&x);

        assert_eq!(dx, Lc::from_iter([
            (KhIGen::Q(
                KhGen::new(
                    State::from([1,1,0]),
                    KhLabel::from([X,X])
                )
            ), R::one()),

            (KhIGen::Q(
                KhGen::new(
                    State::from([0,1,1]),
                    KhLabel::from([X,X])
                )
            ), R::one()),
        ]));
    }

    #[test]
    fn complex_kh() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None).into_complex();

        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 10);
        assert_eq!(c.rank(2), 18);
        assert_eq!(c.rank(3), 20);
        assert_eq!(c.rank(4), 8);

        c.check_d_all();
    }

    #[test]
    fn complex_fbn() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::one();
        let c = KhIComplex::new(&l, &h, None).into_complex();

        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 10);
        assert_eq!(c.rank(2), 18);
        assert_eq!(c.rank(3), 20);
        assert_eq!(c.rank(4), 8);
        
        c.check_d_all();
    }

    #[test]
    fn complex_bn() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        type P = Poly<'H', R>;
        let h = P::variable();

        let c = KhIComplex::new(&l, &h, None).into_complex();

        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 10);
        assert_eq!(c.rank(2), 18);
        assert_eq!(c.rank(3), 20);
        assert_eq!(c.rank(4), 8);
        
        c.check_d_all();
    }

    #[test]
    fn complex_red() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );
        let red_e = Some(3);

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, red_e).into_complex();

        assert_eq!(c.rank(0), 2);
        assert_eq!(c.rank(1), 5);
        assert_eq!(c.rank(2), 9);
        assert_eq!(c.rank(3), 10);
        assert_eq!(c.rank(4), 4);
        
        c.check_d_all();
    }
}