use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use yui::bitseq::{Bit, BitSeq};
use yui::{Ring, RingOps};
use yui_link::State;

use crate::v1::cube::KhCube;
use crate::{KhGen, KhLabel};
use super::InvLink;

pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    cube: KhCube<R>,
    state_map: HashMap<State, State>,
    label_map: HashMap<State, HashMap<usize, usize>>
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R) -> Self { 
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");

        let t = R::zero();
        let cube = KhCube::new(l.link(), h, &t);

        let n = cube.dim();
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

        Self { cube, state_map, label_map }
    }

    fn tau_state(&self, s: State) -> State { 
        self.state_map.get(&s).cloned().unwrap()
    }

    fn tau_label(&self, s: State, l: KhLabel) -> KhLabel { 
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

    fn tau_gen(&self, x: KhGen) -> KhGen {
        let s = x.state;
        let t = self.tau_state(s);
        let l = self.tau_label(s, x.label);
        KhGen::new(t, l)
    }
}

#[cfg(test)]
mod tests {
    use yui::FF;
    use num_traits::Zero;
    use yui_link::Link;
    use super::*;
 
    #[test]
    fn tau_state() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h);

        assert_eq!(
            c.tau_state(State::from([0,0,0])), 
            State::from([0,0,0])
        );
        
        assert_eq!(
            c.tau_state(State::from([1,0,0])), 
            State::from([1,0,0])
        );

        assert_eq!(
            c.tau_state(State::from([0,1,0])), 
            State::from([0,0,1])
        );

        assert_eq!(
            c.tau_state(State::from([0,0,1])), 
            State::from([0,1,0])
        );

        assert_eq!(
            c.tau_state(State::from([1,1,0])), 
            State::from([1,0,1])
        );

        assert_eq!(
            c.tau_state(State::from([1,0,1])), 
            State::from([1,1,0])
        );

        assert_eq!(
            c.tau_state(State::from([0,1,1])), 
            State::from([0,1,1])
        );

        assert_eq!(
            c.tau_state(State::from([1,1,1])), 
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
        let c = KhIComplex::new(&l, &h);

        assert_eq!(
            c.tau_label(State::from([0,1,0]), KhLabel::from([X])), 
            KhLabel::from([X])
        );
    }
}