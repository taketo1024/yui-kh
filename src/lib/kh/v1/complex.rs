use yui::{Ring, RingOps};
use yui_link::Link;

use super::cube::KhCube;
use crate::{KhComplex, KhComplexBigraded};

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v1(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        let red_e = reduced.then(|| l.first_edge().unwrap());
        let deg_shift = Self::deg_shift_for(l, reduced);
        
        let cube = KhCube::new(l, h, t, red_e, deg_shift);
        let complex = cube.into_complex();

        let canon_cycles = if t.is_zero() && l.is_knot() {
            let p = l.first_edge().unwrap();
            Self::make_canon_cycles(l, p, &R::zero(), h, reduced, deg_shift)
        } else { 
            vec![]
        };

        KhComplex::new_impl(complex, canon_cycles, reduced, deg_shift)
    }        
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v1(l: Link, reduced: bool) -> Self { 
        let c = KhComplex::new_v1(&l, &R::zero(), &R::zero(), reduced);
        c.into_bigraded()
    }
}

#[cfg(test)]
mod tests {
    use yui_homology::{RModStr, ChainComplexCommon};
    use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_v1(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);

        assert_eq!(c[-3].rank(), 8);
        assert_eq!(c[-2].rank(), 12);
        assert_eq!(c[-1].rank(), 6);
        assert_eq!(c[ 0].rank(), 4);

        c.check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::trefoil();
        let c = KhComplex::new_v1(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);

        assert_eq!(c[-3].rank(), 4);
        assert_eq!(c[-2].rank(), 6);
        assert_eq!(c[-1].rank(), 3);
        assert_eq!(c[ 0].rank(), 2);

        c.check_d_all();
    }
}
