use std::collections::HashSet;

use yui::bitseq::Bit;
use yui::{Ring, RingOps};
use yui_link::State;

use crate::v1::cube::KhCube;
use super::InvLink;

pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    l: InvLink,
    c0: KhCube<R>, // lower
    c1: KhCube<R>  // upper
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R) -> Self { 
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");

        let l = l.clone();
        let t = R::zero();
        let c0 = KhCube::new(l.link(), h, &t);
        let c1 = KhCube::new(l.link(), h, &t);

        Self { l, c0, c1 }
    }

    fn tau_state(&self, s: State) -> State { 
        let n = s.len();
        let xs = s.iter().enumerate().filter_map(|(i, b)| 
            b.is_one().then_some(i)
        ).map(|i| 
            self.l.inv_x(i)
        ).collect::<HashSet<_>>();

        State::from_iter((0..n).map(|i| 
            if xs.contains(&i) { Bit::Bit1 } else { Bit::Bit0 }
        ))
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
}