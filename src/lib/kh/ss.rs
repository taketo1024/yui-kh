// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use core::panic;

use itertools::Itertools;
use log::info;
use yui_homology::{ChainComplexTrait, RModStr};
use yui_link::Link;
use yui_homology::utils::{ChainReducer, HomologyCalc};
use yui::{EucRing, EucRingOps};

use crate::misc::div_vec;
use crate::kh::KhComplex;

#[derive(Clone, Copy)]
pub enum Ver { 
    V1, V2
}

pub fn ss_invariant<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    ss_invariant_ver(l, c, reduced, Ver::V2)
}

pub fn ss_invariant_v1<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    ss_invariant_ver(l, c, reduced, Ver::V1)
}

pub fn ss_invariant_v2<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    ss_invariant_ver(l, c, reduced, Ver::V2)
}

fn ss_invariant_ver<R>(l: &Link, c: &R, reduced: bool, ver: Ver) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());

    info!("compute ss, c = {c} ({}).", std::any::type_name::<R>());

    let d = match ver {
        Ver::V1 => div_v1(l, c, reduced),
        Ver::V2 => div_v2(l, c, reduced),
    };
    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let ss = 2 * d + w - r + 1;

    info!("d = {d}, w = {w}, r = {r}.");
    info!("ss = {ss} (c = {c}, {}).", if reduced { "reduced" } else { "unreduced" } );

    ss
}

fn div_v1<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let r = if reduced { 1 } else { 2 };

    let ckh = KhComplex::new_v1(l, c, &R::zero(), reduced);
    let zs = ckh.canon_cycles();

    assert_eq!(zs.len(), r);
    assert!(zs.iter().all(|z| z.gens().all(|x| x.h_deg() == 0)));

    let vs = zs.iter().map(|z| {
        ckh[0].vectorize(z)
    }).collect_vec();

    let range = -1 ..= 1;
    let mut reducer = ChainReducer::new(range.clone(), 1);

    for i in range { 
        reducer.set_matrix(i, ckh.d_matrix(i), false);
    }

    for v in vs { 
        reducer.add_vec(0, v);
    }

    reducer.reduce_all(false);
    reducer.reduce_all(true);

    let dm = reducer.matrix(-1).unwrap().clone();
    let d0 = reducer.matrix( 0).unwrap().clone();
    let vs = reducer.vecs(0).unwrap().clone();

    let kh = HomologyCalc::calculate(dm, d0, true);

    info!("homology: {}", kh.math_symbol());    
    assert_eq!(kh.rank(), r);

    let t = kh.trans().unwrap();
    let vs = vs.iter().map(|v| {
        t.forward(v).subvec(0..r)
    }).collect_vec();

    let v = &vs[0];
    let Some(d) = div_vec(v, c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v, c)
    };
    
    d
}

fn div_v2<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let r = if reduced { 1 } else { 2 };

    let ckh = KhComplex::new_v2(l, c, &R::zero(), reduced);
    let d0 = ckh.d_matrix(-1);
    let d1 = ckh.d_matrix( 0);
    let zs = ckh.canon_cycles();

    assert_eq!(zs.len(), r);
    assert!(zs.iter().all(|z| z.gens().all(|x| x.h_deg() == 0)));

    let kh = HomologyCalc::calculate(d0, d1, true);

    info!("homology: {}", kh.math_symbol());    
    assert_eq!(kh.rank(), r);
    
    let t = kh.trans().unwrap();
    let vs = zs.into_iter().map(|z| {
        let v = ckh[0].vectorize(&z);
        let v = t.forward(&v).subvec(0..r);
        info!("{z} -> {v}");
        v
    }).collect_vec();

    let v = &vs[0];
    let Some(d) = div_vec(v, c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v, c)
    };
    
    d
}

#[cfg(test)]
mod tests {
    use yui_link::Link;
    use super::*;

    fn test(l: &Link, c: i32, value: i32) { 
        for r in [true, false] { 
            assert_eq!(ss_invariant(l, &c, r), value);
            assert_eq!(ss_invariant(&l.mirror(), &c, r), -value);
        }
    }

    fn test_v1(l: &Link, c: i32, value: i32) { 
        for r in [true, false] { 
            assert_eq!(ss_invariant_v1(l, &c, r), value);
            assert_eq!(ss_invariant_v1(&l.mirror(), &c, r), -value);
        }
    }

    #[test]
    fn test_unknot() { 
        let l = Link::unknot();
        test(&l, 2, 0);
    }

    #[test]
    fn test_unknot_v1() { 
        let l = Link::unknot();
        test_v1(&l, 2, 0);
    }

    #[test]
    fn test_unknot_rm1() { 
        let l = Link::from_pd_code([[0,0,1,1]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_unknot_rm1_neg() { 
        let l = Link::from_pd_code([[0,1,1,0]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_3_1() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_4_1() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_5_1() { 
        let l = Link::from_pd_code([[1,6,2,7],[3,8,4,9],[5,10,6,1],[7,2,8,3],[9,4,10,5]]);
        test(&l, 2, -4);
    }

    #[test]
    fn test_5_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,8,4,9],[5,10,6,1],[9,6,10,7],[7,2,8,3]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_6_1() { 
        let l = Link::from_pd_code([[1,4,2,5],[7,10,8,11],[3,9,4,8],[9,3,10,2],[5,12,6,1],[11,6,12,7]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_6_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[5,10,6,11],[3,9,4,8],[9,3,10,2],[7,12,8,1],[11,6,12,7]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_6_3() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_7_1() { 
        let l = Link::from_pd_code([[1,8,2,9],[3,10,4,11],[5,12,6,13],[7,14,8,1],[9,2,10,3],[11,4,12,5],[13,6,14,7]]);
        test(&l, 2, -6);
    }

    #[test]
    fn test_7_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,10,4,11],[5,14,6,1],[7,12,8,13],[11,8,12,9],[13,6,14,7],[9,2,10,3]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_7_3() { 
        let l = Link::from_pd_code([[6,2,7,1],[10,4,11,3],[14,8,1,7],[8,14,9,13],[12,6,13,5],[2,10,3,9],[4,12,5,11]]);
        test(&l, 2, 4);
    }

    #[test]
    fn test_8_19() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        test(&l, 2, 6);
    }

    #[test]
    fn test_8_19_v1() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        test_v1(&l, 2, 6);
    }

    #[test]
    #[ignore]
    fn test_k14() { 
        let l = Link::from_pd_code([[1,19,2,18],[19,1,20,28],[20,13,21,14],[12,17,13,18],[16,21,17,22],[5,15,6,14],[15,5,16,4],[6,27,7,28],[2,7,3,8],[26,3,27,4],[25,23,26,22],[11,9,12,8],[23,10,24,11],[9,24,10,25]]);
        assert_eq!(ss_invariant_v1(&l, &2, true), -2);
        assert_eq!(ss_invariant_v1(&l, &3, true), 0);
    }
}