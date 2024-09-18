use yui::{Ring, RingOps};
use yui_link::Link;

use crate::kh::{KhComplex, KhComplexBigraded};

use super::builder::TngComplexBuilder;

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v2(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        TngComplexBuilder::build_kh_complex(l, h, t, reduced)
    }        
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v2(l: Link, reduced: bool) -> Self { 
        let c = KhComplex::new_v2(&l, &R::zero(), &R::zero(), reduced);
        c.into_bigraded()
    }
}

#[cfg(test)]
mod tests {
    use yui_homology::{ChainComplexCommon, RModStr};
    use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_v2(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);

        assert_eq!(c[-3].rank(), 2);
        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);

        c.check_d_all();
    }
}
