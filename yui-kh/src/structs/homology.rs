use std::ops::{RangeInclusive, Index};
use delegate::delegate;

use yui_homology::{GridTrait, XHomology, XHomologySummand, XHomology2};
use yui::{EucRing, EucRingOps, isize2};
use yui_link::Link;

use crate::KhEnhState;

pub struct KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: XHomology<KhEnhState, R>
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        KhHomology::new_v2(l, h, t, reduced)
    }

    pub(crate) fn _new(inner: XHomology<KhEnhState, R>) -> Self { 
        Self { inner }
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        let h = &self.inner;
        let h_min = h.support().min().unwrap_or(0);
        let h_max = h.support().max().unwrap_or(0);
        h_min ..= h_max
    }
}

impl<R> GridTrait<isize> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Itr = std::vec::IntoIter<isize>;
    type E = XHomologySummand<KhEnhState, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::E;
        }
    }
}

impl<R> Index<isize> for KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = XHomologySummand<KhEnhState, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

pub struct KhHomologyBigraded<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: XHomology2<KhEnhState, R>
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: Link, reduced: bool) -> Self {
        KhHomologyBigraded::new_v2(l, reduced)
    }

    pub(crate) fn _new(inner: XHomology2<KhEnhState, R>) -> Self { 
        Self { inner }
    }
}

impl<R> GridTrait<isize2> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Itr = std::vec::IntoIter<isize2>;
    type E = XHomologySummand<KhEnhState, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize2) -> bool;
            fn get(&self, i: isize2) -> &Self::E;
        }
    }
}

impl<R> Index<(isize, isize)> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = XHomologySummand<KhEnhState, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

#[cfg(test)]
mod tests {
    use yui_homology::RModStr;
    use yui_link::Link;
    use super::*;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -3..=0);

        assert_eq!(h[-3].rank(), 1);
        assert_eq!(h[-3].is_free(), true);

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].tors(), &vec![2]);

        assert_eq!(h[-1].is_zero(), true);

        assert_eq!(h[ 0].rank(), 2);
        assert_eq!(h[ 0].is_free(), true);
    }

    #[test]
    fn kh_trefoil_mirror() {
        let l = Link::trefoil().mirror();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=3);

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].is_zero(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -2..=2);

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].is_free(), true);

        assert_eq!(h[-1].rank(), 1);
        assert_eq!(h[-1].tors(), &vec![2]);

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 1);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].tors(), &vec![2]);
    }

    #[test]
    fn kh_empty_bigr() {
        let l = Link::empty();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(0,0)].rank(), 1);
        assert_eq!(h[(0,0)].is_free(), true);
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(0,-1)].rank(), 1);
        assert_eq!(h[(0,-1)].is_free(), true);
        assert_eq!(h[(0, 1)].rank(), 1);
        assert_eq!(h[(0, 1)].is_free(), true);
    }

    #[test]
    fn kh_unknot_bigr_red() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h[(0, 0)].is_free(), true);
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(-3,-9)].rank(), 1);
        assert_eq!(h[(-3,-9)].is_free(), true);
        assert_eq!(h[(-2,-7)].rank(), 0);
        assert_eq!(h[(-2,-7)].tors(), &vec![2]);
        assert_eq!(h[(-2,-5)].rank(), 1);
        assert_eq!(h[(-2,-5)].is_free(), true);
        assert_eq!(h[( 0,-3)].rank(), 1);
        assert_eq!(h[( 0,-3)].is_free(), true);
        assert_eq!(h[( 0,-1)].rank(), 1);
        assert_eq!(h[( 0,-1)].is_free(), true);
    }

    #[test]
    fn kh_trefoil_mirror_bigr() {
        let l = Link::trefoil().mirror();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(0, 1)].rank(), 1);
        assert_eq!(h[(0, 1)].is_free(), true);
        assert_eq!(h[(0, 3)].rank(), 1);
        assert_eq!(h[(0, 3)].is_free(), true);
        assert_eq!(h[(2, 5)].rank(), 1);
        assert_eq!(h[(2, 5)].is_free(), true);
        assert_eq!(h[(3, 7)].rank(), 0);
        assert_eq!(h[(3, 7)].tors(), &vec![2]);
        assert_eq!(h[(3, 9)].rank(), 1);
        assert_eq!(h[(3, 9)].is_free(), true);
    }

    #[test]
    fn kh_trefoil_bigr_red() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[(-3,-8)].rank(), 1);
        assert_eq!(h[(-3,-8)].is_free(), true);
        assert_eq!(h[(-2,-6)].rank(), 1);
        assert_eq!(h[(-2,-6)].is_free(), true);
        assert_eq!(h[( 0,-2)].rank(), 1);
        assert_eq!(h[( 0,-2)].is_free(), true);
    }

    #[test]
    fn kh_figure8_bigr() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(-2,-5)].rank(), 1);
        assert_eq!(h[(-2,-5)].is_free(), true);
        assert_eq!(h[(-1,-3)].rank(), 0);
        assert_eq!(h[(-1,-3)].tors(), &vec![2]);
        assert_eq!(h[(-1,-1)].rank(), 1);
        assert_eq!(h[(-1,-1)].is_free(), true);
        assert_eq!(h[( 0,-1)].rank(), 1);
        assert_eq!(h[( 0,-1)].is_free(), true);
        assert_eq!(h[( 0, 1)].rank(), 1);
        assert_eq!(h[( 0, 1)].is_free(), true);
        assert_eq!(h[( 1, 1)].rank(), 1);
        assert_eq!(h[( 1, 1)].is_free(), true);
        assert_eq!(h[( 2, 3)].rank(), 0);
        assert_eq!(h[( 2, 3)].tors(), &vec![2]);
        assert_eq!(h[( 2, 5)].rank(), 1);
        assert_eq!(h[( 2, 5)].is_free(), true);
   }

    #[test]
    fn kh_figure8_bigr_red() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[(-2,-4)].rank(), 1);
        assert_eq!(h[(-2,-4)].is_free(), true);
        assert_eq!(h[(-1,-2)].rank(), 1);
        assert_eq!(h[(-1,-2)].is_free(), true);
        assert_eq!(h[( 0, 0)].rank(), 1);
        assert_eq!(h[( 0, 0)].is_free(), true);
        assert_eq!(h[( 1, 2)].rank(), 1);
        assert_eq!(h[( 1, 2)].is_free(), true);
        assert_eq!(h[( 2, 4)].rank(), 1);
        assert_eq!(h[( 2, 4)].is_free(), true);
   }
}