// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use core::panic;

use log::info;
use yui_homology::RModStr;
use yui::{EucRing, EucRingOps};
use yui_link::InvLink;

use crate::misc::div_vec;
use crate::khi::KhIComplex;

#[allow(unused)]
pub fn ssi_invariants<R>(l: &InvLink, c: &R) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());
    assert!(l.link().is_knot());

    info!("compute ssi, c = {c} over {}.", R::math_symbol());

    let w = l.link().writhe();
    let r = l.link().seifert_circles().len() as i32;
    let (d0, d1) = div(l, c);

    let ss0 = 2 * d0 + w - r + 1;
    let ss1 = 2 * d1 + w - r + 1;

    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}.");
    info!("ssi = ({ss0}, {ss1}).");

    (ss0, ss1)
}

fn div<R>(l: &InvLink, c: &R) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let ckh = KhIComplex::new(l, c, true).truncated(-1..=2);
    let zs = ckh.canon_cycles();

    assert_eq!(zs.len(), 2);

    let z0 = &zs[0];
    let z1 = &zs[1]; 

    assert!(z0.gens().all(|x| x.h_deg() == 0));
    assert!(z1.gens().all(|x| x.h_deg() == 1));

    let kh = ckh.homology();
    assert_eq!(kh[0].rank(), 1);
    assert_eq!(kh[1].rank(), 1);    

    info!("KhI[0]: {}", kh[0].math_symbol());    
    info!("KhI[1]: {}", kh[1].math_symbol());    

    let v0 = kh[0].vectorize(&z0).subvec(0..1);
    let v1 = kh[1].vectorize(&z1).subvec(0..1);

    info!("a0: {z0} -> [{}]", v0.to_dense()[0]);
    info!("a1: {z1} -> [{}]", v1.to_dense()[0]);

    let Some(d0) = div_vec(&v0, c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v0, c)
    };
    let Some(d1) = div_vec(&v1, c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v1, c)
    };

    (d0, d1)
}

// TODO 
#[cfg(all(test, feature = "old"))]
mod tests {
    use yui::poly::HPoly;
    use yui::FF2;
    use yui_link::Link;

    use super::*;

    type R = FF2;
    type P = HPoly<'H', R>;

    #[test]
    fn test_unknot_pos_twist() { 
        let l = InvLink::new(
            Link::from_pd_code([[0,0,1,1]]),
            [],
            Some(0)
        );
        let c = P::variable();

        let div = div(&l, &c);
        assert_eq!(div.0, 0);
        assert_eq!(div.1, 0);

        let ssi = ssi_invariants(&l, &c);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist() { 
        let l = InvLink::new(
            Link::from_pd_code([[0,1,1,0]]),
            [],
            Some(0)
        );
        let c = P::variable();

        let div = div(&l, &c);
        assert_eq!(div.0, 1);
        assert_eq!(div.1, 1);

        let ssi = ssi_invariants(&l, &c);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist2() { 
        let l = InvLink::new(
            Link::from_pd_code([[0,1,3,0],[2,3,1,2]]),
            [(1,3)],
            Some(0)
        );
        let c = P::variable();

        let div = div(&l, &c);
        assert_eq!(div.0, 2);
        assert_eq!(div.1, 2);

        let ssi = ssi_invariants(&l, &c);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_3_1() { 
        let l = InvLink::new( // positive
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            Some(3)
        );
        let c = P::variable();

        let div = div(&l, &c);
        assert_eq!(div.0, 0);
        assert_eq!(div.1, 0);

        let ssi = ssi_invariants(&l, &c);
        assert_eq!(ssi.0, 2);
        assert_eq!(ssi.1, 2);
    }

    #[test]
    fn test_3_1_m() { 
        let l = InvLink::new( // negative
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            Some(3)
        ).mirror();
        let c = P::variable();

        let div = div(&l, &c);
        assert_eq!(div.0, 1);
        assert_eq!(div.1, 1);

        let ssi = ssi_invariants(&l, &c);
        assert_eq!(ssi.0, -2);
        assert_eq!(ssi.1, -2);
    }
}