use std::collections::{HashMap, HashSet};
use std::ops::RangeInclusive;

use itertools::Itertools;
use yui::bitseq::{Bit, BitSeq};
use yui::lc::Lc;
use yui::{Ring, RingOps};
use yui_homology::{Grid, XChainComplex, XModStr};
use yui_link::{InvLink, State};

use crate::v1::cube::KhCube;
use crate::{KhGen, KhLabel};
use super::KhIGen;

pub struct KhICube<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    cube: KhCube<R>,
    state_map: HashMap<State, State>,
    label_map: HashMap<State, HashMap<usize, usize>>,
    deg_shift: (isize, isize)
}

impl<R> KhICube<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R, reduced: bool, deg_shift: (isize, isize)) -> Self { 
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        assert!(!reduced || l.base_pt().is_some());

        let n = l.link().crossing_num();
        let t = R::zero();
        let reduce_e = if reduced { l.base_pt() } else { None };
        let cube = KhCube::new(l.link(), h, &t, reduce_e, deg_shift);

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
        KhGen::new(t, l, x.deg_shift)
    }

    // f = 1 + τ
    fn f(&self, x: &KhGen) -> Lc<KhGen, R> { 
        let x = x.clone();
        let tx = self.t(&x);
        Lc::from_iter([
            (x,  R::one()), 
            (tx, R::one())
        ])
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        *self.cube.h_range().start() ..= self.cube.h_range().end() + 1
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        self.cube.q_range()
    }

    pub fn generators(&self, i: isize) -> Vec<KhIGen> { 
        Iterator::chain(
            self.cube.generators(i).iter().map(|&&x| 
                KhIGen::B(x)
            ),
            self.cube.generators(i - 1).iter().map(|&&x|
                KhIGen::Q(x)
            )
        ).collect()
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
                let dx = self.cube.d(x).map_gens(|&y| 
                    KhIGen::B(y)
                );
                let fx = self.f(x).map_gens(|&y| 
                    KhIGen::Q(y)
                );
                dx + fx
            },
            KhIGen::Q(x) => {
                self.cube.d(x).map_gens(|&y| 
                    KhIGen::Q(y)
                )
            },
        }
    }

    pub fn into_complex(self) -> XChainComplex<KhIGen, R> {
        XChainComplex::new(
            Grid::generate(self.h_range(), |i| self.summand(i)),
            1, 
            move |_, z| self.differentiate(z)
        )
    }   
}

#[cfg(test)]
mod tests {
    #![allow(unused)]

    use yui::poly::Poly;
    use yui::FF2;
    use num_traits::{Zero, One};
    use yui_homology::{ChainComplexCommon, ChainComplexTrait, DisplaySeq};
    use yui_link::Link;
    use crate::KhHomology;

    use super::*;
 
    #[test]
    fn tau_state() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhICube::new(&l, &h, false, (0, 0));

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
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhICube::new(&l, &h, false, (0, 0));

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
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhICube::new(&l, &h, false, (0, 0));

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
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhICube::new(&l, &h, false, (0, 0));

        let x = KhIGen::B(
            KhGen::new(
                State::from([0,0,0]),
                KhLabel::from([X, I]),
                c.deg_shift
            )
        );
        let dx = c.d(&x);

        assert_eq!(dx, Lc::from_iter([
            (KhIGen::B(
                KhGen::new(
                    State::from([1,0,0]),
                    KhLabel::from([X]),
                    c.deg_shift
                )
            ), R::one()),

            (KhIGen::B(
                KhGen::new(
                    State::from([0,1,0]),
                    KhLabel::from([X]),
                    c.deg_shift
                )
            ), R::one()),
            
            (KhIGen::B(
                KhGen::new(
                    State::from([0,0,1]),
                    KhLabel::from([X]),
                    c.deg_shift
                )
            ), R::one())
        ]));
    }

    #[test]
    fn d_lower_asym() { 
        use crate::KhAlgGen::*;

        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhICube::new(&l, &h, false, (0, 0));

        let x = KhIGen::B(
            KhGen::new(
                State::from([0,1,0]),
                KhLabel::from([X]),
                c.deg_shift
            )
        );
        let dx = c.d(&x);

        assert_eq!(dx, Lc::from_iter([
            (KhIGen::B(
                KhGen::new(
                    State::from([1,1,0]),
                    KhLabel::from([X,X]),
                    c.deg_shift
                )
            ), R::one()),

            (KhIGen::B(
                KhGen::new(
                    State::from([0,1,1]),
                    KhLabel::from([X,X]),
                    c.deg_shift
                )
            ), R::one()),
            
            (KhIGen::Q(
                KhGen::new(
                    State::from([0,1,0]),
                    KhLabel::from([X]),
                    c.deg_shift
                )
            ), R::one()),

                        
            (KhIGen::Q(
                KhGen::new(
                    State::from([0,0,1]),
                    KhLabel::from([X]),
                    c.deg_shift
                )
            ), R::one())
        ]));
    
    }

    #[test]
    fn d_upper() { 
        use crate::KhAlgGen::*;

        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhICube::new(&l, &h, false, (0, 0));

        let x = KhIGen::Q(
            KhGen::new(
                State::from([0,1,0]),
                KhLabel::from([X]),
                c.deg_shift
            )
        );
        let dx = c.d(&x);

        assert_eq!(dx, Lc::from_iter([
            (KhIGen::Q(
                KhGen::new(
                    State::from([1,1,0]),
                    KhLabel::from([X,X]),
                    c.deg_shift
                )
            ), R::one()),

            (KhIGen::Q(
                KhGen::new(
                    State::from([0,1,1]),
                    KhLabel::from([X,X]),
                    c.deg_shift
                )
            ), R::one()),
        ]));
    }
}