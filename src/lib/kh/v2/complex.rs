use itertools::Itertools;
use num_traits::Zero;
use yui::{Ring, RingOps};
use yui_homology::ChainComplexTrait;
use yui_link::Link;

use crate::{KhComplex, KhComplexBigraded};

use super::builder::TngComplexBuilder;

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v2(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        let deg_shift = Self::deg_shift_for(l, reduced);
        let mut b = TngComplexBuilder::new(l, h, t, reduced);

        if t.is_zero() && l.is_knot() {
            b.make_canon_cycles();
        }
        
        b.process();

        let canon_cycles = b.canon_cycles().iter().map(|z| 
            z.eval(h, t, deg_shift)
        ).collect_vec();
        let complex = b.into_complex().eval(h, t);

        assert!(canon_cycles.iter().all(|z| 
            complex.d(0, z).is_zero()
        ));

        Self::new_impl(complex, canon_cycles, reduced, deg_shift)
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
