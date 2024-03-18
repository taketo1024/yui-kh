use yui::{EucRing, EucRingOps};
use yui_link::Link;

use crate::{KhComplex, KhHomology, KhComplexBigraded, KhHomologyBigraded};

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new_v1(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        KhComplex::new_v1(l, h, t, reduced).homology(false)
    }
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new_v1(l: Link, reduced: bool) -> Self {
        KhComplexBigraded::new_v1(l, reduced).homology(false)
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
        let h = KhHomology::new_v1(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::new_v1(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::new_v1(&l, &0, &0, false);

        assert_eq!(h.h_range(), -3..=0);

        assert_eq!(h[-3].rank(), 1);
        assert!(h[-3].is_free());

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].tors(), &vec![2]);

        assert!(h[-1].is_zero());

        assert_eq!(h[ 0].rank(), 2);
        assert!(h[ 0].is_free());
    }

    #[test]
    fn kh_trefoil_mirror() {
        let l = Link::trefoil().mirror();
        let h = KhHomology::new_v1(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=3);

        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());

        assert!(h[1].is_zero());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let h = KhHomology::new_v1(&l, &0, &0, false);

        assert_eq!(h.h_range(), -2..=2);

        assert_eq!(h[-2].rank(), 1);
        assert!(h[-2].is_free());

        assert_eq!(h[-1].rank(), 1);
        assert_eq!(h[-1].tors(), &vec![2]);

        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].tors(), &vec![2]);
    }

    #[test]
    fn kh_empty_bigr() {
        let l = Link::empty();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(0,0)].rank(), 1);
        assert!(h[(0,0)].is_free());
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(0,-1)].rank(), 1);
        assert!(h[(0,-1)].is_free());
        assert_eq!(h[(0, 1)].rank(), 1);
        assert!(h[(0, 1)].is_free());
    }

    #[test]
    fn kh_unknot_bigr_red() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[(0, 0)].rank(), 1);
        assert!(h[(0, 0)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(-3,-9)].rank(), 1);
        assert!(h[(-3,-9)].is_free());
        assert_eq!(h[(-2,-7)].rank(), 0);
        assert_eq!(h[(-2,-7)].tors(), &vec![2]);
        assert_eq!(h[(-2,-5)].rank(), 1);
        assert!(h[(-2,-5)].is_free());
        assert_eq!(h[( 0,-3)].rank(), 1);
        assert!(h[( 0,-3)].is_free());
        assert_eq!(h[( 0,-1)].rank(), 1);
        assert!(h[( 0,-1)].is_free());
    }

    #[test]
    fn kh_trefoil_mirror_bigr() {
        let l = Link::trefoil().mirror();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(0, 1)].rank(), 1);
        assert!(h[(0, 1)].is_free());
        assert_eq!(h[(0, 3)].rank(), 1);
        assert!(h[(0, 3)].is_free());
        assert_eq!(h[(2, 5)].rank(), 1);
        assert!(h[(2, 5)].is_free());
        assert_eq!(h[(3, 7)].rank(), 0);
        assert_eq!(h[(3, 7)].tors(), &vec![2]);
        assert_eq!(h[(3, 9)].rank(), 1);
        assert!(h[(3, 9)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr_red() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[(-3,-8)].rank(), 1);
        assert!(h[(-3,-8)].is_free());
        assert_eq!(h[(-2,-6)].rank(), 1);
        assert!(h[(-2,-6)].is_free());
        assert_eq!(h[( 0,-2)].rank(), 1);
        assert!(h[( 0,-2)].is_free());
    }

    #[test]
    fn kh_figure8_bigr() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[(-2,-5)].rank(), 1);
        assert!(h[(-2,-5)].is_free());
        assert_eq!(h[(-1,-3)].rank(), 0);
        assert_eq!(h[(-1,-3)].tors(), &vec![2]);
        assert_eq!(h[(-1,-1)].rank(), 1);
        assert!(h[(-1,-1)].is_free());
        assert_eq!(h[( 0,-1)].rank(), 1);
        assert!(h[( 0,-1)].is_free());
        assert_eq!(h[( 0, 1)].rank(), 1);
        assert!(h[( 0, 1)].is_free());
        assert_eq!(h[( 1, 1)].rank(), 1);
        assert!(h[( 1, 1)].is_free());
        assert_eq!(h[( 2, 3)].rank(), 0);
        assert_eq!(h[( 2, 3)].tors(), &vec![2]);
        assert_eq!(h[( 2, 5)].rank(), 1);
        assert!(h[( 2, 5)].is_free());
   }

    #[test]
    fn kh_figure8_bigr_red() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[(-2,-4)].rank(), 1);
        assert!(h[(-2,-4)].is_free());
        assert_eq!(h[(-1,-2)].rank(), 1);
        assert!(h[(-1,-2)].is_free());
        assert_eq!(h[( 0, 0)].rank(), 1);
        assert!(h[( 0, 0)].is_free());
        assert_eq!(h[( 1, 2)].rank(), 1);
        assert!(h[( 1, 2)].is_free());
        assert_eq!(h[( 2, 4)].rank(), 1);
        assert!(h[( 2, 4)].is_free());
   }
}