use std::ops::Index;
use delegate::delegate;

use yui::lc::Lc;
use yui::{EucRing, EucRingOps, Ring, RingOps};
use yui_homology::{ChainComplexTrait, GridTrait, XChainComplex, XChainComplexSummand, XHomology};
use yui_link::Edge;
use yui_matrix::sparse::SpMat;

use crate::KhComplex;

use super::{InvLink, KhIGen, KhICube};

pub type KhIChain<R> = Lc<KhIGen, R>;
pub type KhIComplexSummand<R> = XChainComplexSummand<KhIGen, R>;

pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    inner: XChainComplex<KhIGen, R>,
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R, reduce_e: Option<Edge>) -> Self { 
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduce_e.is_some());
        let cube = KhICube::new(l, h, reduce_e, deg_shift);
        let inner = cube.into_complex();

        Self { inner }
    }
}

impl<R> Index<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhIComplexSummand<R>;

    delegate! { 
        to self.inner {
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<isize>;
    type Output = KhIComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Output;
        }
    }
}

impl<R> ChainComplexTrait<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Element = KhIChain<R>;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize) -> usize;
            fn d_deg(&self) -> isize;
            fn d(&self, i: isize, z: &Self::Element) -> Self::Element;
            fn d_matrix(&self, i: isize) -> SpMat<R>;
        }
    }
}

impl<R> KhIComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(self, with_trans: bool) -> XHomology<KhIGen, R> {
        self.inner.reduced().homology(with_trans)
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
    fn complex_kh() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)]
        );

        type R = FF<2>;
        let h = R::zero();
        let c = KhIComplex::new(&l, &h, None);

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
        let c = KhIComplex::new(&l, &h, None);

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

        let c = KhIComplex::new(&l, &h, None);

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
        let c = KhIComplex::new(&l, &h, red_e);

        assert_eq!(c.rank(0), 2);
        assert_eq!(c.rank(1), 5);
        assert_eq!(c.rank(2), 9);
        assert_eq!(c.rank(3), 10);
        assert_eq!(c.rank(4), 4);
        
        c.check_d_all();
    }
}