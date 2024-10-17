// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use core::panic;
use itertools::Itertools;
use num_traits::Zero;
use log::info;

use yui_homology::RModStr;
use yui::{EucRing, EucRingOps};
use yui_link::InvLink;

use crate::kh::KhChainExt;
use crate::misc::div_vec;
use crate::khi::KhIComplex;

pub fn ssi_invariants<R>(l: &InvLink, c: &R, reduced: bool) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());
    assert!(l.link().is_knot());

    info!("compute ssi, c = {c} over {}.", R::math_symbol());

    let w = l.link().writhe();
    let r = l.link().seifert_circles().len() as i32;
    let (d0, d1) = div(l, c, reduced);

    let ss0 = 2 * d0 + w - r + 1;
    let ss1 = 2 * d1 + w - r + 1;

    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}.");
    info!("ssi = ({ss0}, {ss1}).");

    (ss0, ss1)
}

fn div<R>(l: &InvLink, c: &R, reduced: bool) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let r = if reduced { 1 } else { 2 };
    let t = R::zero(); 
    let ckh = KhIComplex::new(l, c, &t, reduced).truncated(-1..=2);
    let zs = ckh.canon_cycles();

    assert_eq!(zs.len(), 2 * r);

    for (i, z) in zs.iter().enumerate() { 
        let h = z.h_deg();
        let v = ckh[h].vectorize(z);
        info!("a[{i}] in CKh[{h}]: {z} -> [{}]", v.to_dense().iter().join(","));
    }

    assert!(zs.iter().all(|z| !z.is_zero()));

    let z0 = &zs[0];
    let z1 = &zs[r]; 

    assert!(z0.gens().all(|x| x.h_deg() == 0));
    assert!(z1.gens().all(|x| x.h_deg() == 1));

    let kh = ckh.homology();
    assert_eq!(kh[0].rank(), r);
    assert_eq!(kh[1].rank(), r);    

    info!("KhI[0]: {}", kh[0].math_symbol());    
    info!("KhI[1]: {}", kh[1].math_symbol());    

    let v0 = kh[0].vectorize(&z0).subvec(0..r);
    let v1 = kh[1].vectorize(&z1).subvec(0..r);

    info!("a0: {z0} -> [{}]", v0.to_dense().iter().join(","));
    info!("a1: {z1} -> [{}]", v1.to_dense().iter().join(","));

    let Some(d0) = div_vec(&v0, c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v0, c)
    };
    let Some(d1) = div_vec(&v1, c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v1, c)
    };

    (d0, d1)
}

#[cfg(test)]
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

        let ssi = ssi_invariants(&l, &c, false);
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

        let ssi = ssi_invariants(&l, &c, false);
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

        let ssi = ssi_invariants(&l, &c, true);
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

        let ssi = ssi_invariants(&l, &c, true);
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

        let ssi = ssi_invariants(&l, &c, true);
        assert_eq!(ssi.0, -2);
        assert_eq!(ssi.1, -2);
    }

    macro_rules! test {
        ($(#[$m:meta])* $test:ident, $name:literal, $expected:expr) => {
            $(#[$m])* 
            #[test]
            fn $test() -> Result<(), Box<dyn std::error::Error>> { 
                type R = FF2;
                type P = HPoly<'H', R>;
                let c = P::variable();
    
                let l = InvLink::load($name)?;
                let ssi = ssi_invariants(&l, &c, true);
                assert_eq!(ssi, $expected);
    
                Ok(())
            }
        }
    }
    
    test!(k3_1, "3_1", (2, 2));
    test!(k4_1, "4_1", (0, 0));
    test!(k5_1, "5_1", (4, 4));
    test!(k5_2a, "5_2a", (2, 2));
    test!(k5_2b, "5_2b", (2, 2));
    test!(k6_1a, "6_1a", (0, 0));
    test!(k6_1b, "6_1b", (0, 0));
    test!(k6_2a, "6_2a", (2, 2));
    test!(k6_2b, "6_2b", (2, 2));
    test!(k6_3, "6_3", (0, 0));
    test!(k7_1, "7_1", (6, 6));
    test!(k7_2a, "7_2a", (2, 2));
    test!(k7_2b, "7_2b", (2, 2));
    test!(k7_3a, "7_3a", (4, 4));
    test!(k7_3b, "7_3b", (4, 4));
    test!(k7_4a, "7_4a", (2, 2));
    test!(k7_4b, "7_4b", (2, 2));
    test!(k7_5a, "7_5a", (4, 4));
    test!(k7_5b, "7_5b", (4, 4));
    test!(k7_6a, "7_6a", (-2, -2));
    test!(k7_6b, "7_6b", (-2, -2));
    test!(k7_7a, "7_7a", (0, 0));
    test!(k7_7b, "7_7b", (0, 0));

    #[test]
    fn k9_46() { 
        let l = InvLink::sinv_knot_from_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );

        let c = P::variable();
        let ssi = ssi_invariants(&l, &c, true);

        assert_eq!(ssi, (0, 2));
    }   
}