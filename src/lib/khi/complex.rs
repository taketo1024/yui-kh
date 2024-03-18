use yui::{Ring, RingOps};

use crate::v1::cube::KhCube;
use super::InvLink;

pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    c0: KhCube<R>, // lower
    c1: KhCube<R>  // upper
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R) -> Self { 
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        let t = R::zero();
        let c0 = KhCube::new(l.link(), h, &t);
        let c1 = KhCube::new(l.link(), h, &t);
        Self { c0, c1 }
    }
}

#[cfg(test)]
mod tests {
    use yui::FF;
    use num_traits::Zero;
    use yui_link::Link;
    use super::*;
 
    #[test]
    fn init() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h);
    }
}