use std::ops::{Index, RangeInclusive};
use cartesian::cartesian;
use delegate::delegate;

use yui::lc::Lc;
use yui::{EucRing, EucRingOps, Ring, RingOps};
use yui_homology::{isize2, ChainComplexTrait, Grid1, Grid2, GridTrait, XChainComplex, XChainComplex2, XChainComplexSummand, XModStr};
use yui_link::InvLink;
use yui_matrix::sparse::SpMat;

use crate::kh::{KhChain, KhComplex, KhGen};
use crate::khi::KhIHomology;
use crate::khi::KhIGen;

pub type KhIChain<R> = Lc<KhIGen, R>;

pub trait KhIChainExt { 
    fn h_deg(&self) -> isize;
    fn q_deg(&self) -> isize;
}

impl<R> KhIChainExt for KhIChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn h_deg(&self) -> isize {
        self.gens().map(|x| x.h_deg()).min().unwrap_or(0)
    }
    
    fn q_deg(&self) -> isize {
        self.gens().map(|x| x.q_deg()).min().unwrap_or(0)
    }
}

pub type KhIComplexSummand<R> = XChainComplexSummand<KhIGen, R>;

pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    inner: XChainComplex<KhIGen, R>,
    canon_cycles: Vec<KhIChain<R>>,
    deg_shift: (isize, isize)
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R, reduced: bool) -> Self { 
        Self::new_v2(l, h, reduced)
    }

    pub fn new_v2(l: &InvLink, h: &R, reduced: bool) -> Self {
        use crate::khi::internal::v2::builder::SymTngBuilder;

        let t = R::zero(); // TODO
        SymTngBuilder::build_khi_complex(l, h, &t, reduced)
    }

    #[cfg(feature = "old")]
    pub fn new_v1(l: &InvLink, h: &R, reduced: bool) -> Self { 
        use crate::khi::v1::cube::KhICube;
        use crate::kh::KhComplex;

        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        assert!(!reduced || l.base_pt().is_some());

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduced);
        let cube = KhICube::new(l, h, reduced, deg_shift);
        let inner = cube.into_complex();

        let canon_cycles = if l.base_pt().is_some() && l.link().is_knot() {
            let p = l.base_pt().unwrap();
            let zs = KhComplex::make_canon_cycles(l.link(), p, &R::zero(), h, reduced, deg_shift);
            Iterator::chain(
                zs.iter().map(|z| z.map_gens(|x| KhIGen::B(*x))),
                zs.iter().map(|z| z.map_gens(|x| KhIGen::Q(*x)))
            ).collect()
        } else { 
            vec![]
        };

        Self::new_impl(inner, canon_cycles, deg_shift)
    }

    pub fn from_kh_complex<'a, F>(c: KhComplex<R>, map: F) -> Self
    where F: Fn(&KhGen) -> KhGen + Send + Sync + 'static {
        let deg_shift = c.deg_shift();
        let h_range = c.h_range();
        let h_range = *h_range.start() ..= (h_range.end() + 1);

        let summands = Grid1::generate(h_range, |i| { 
            let b_gens = c[i].gens().iter().map(|x| KhIGen::B(*x));
            let q_gens = c[i - 1].gens().iter().map(|x| KhIGen::Q(*x));
            XModStr::free(b_gens.chain(q_gens))
        });

        let d = move |i: isize, x: &KhIGen| -> KhIChain<R> { 
            match x { 
                KhIGen::B(x) => {
                    let z = KhChain::from(*x);
                    let dx = c.d(i, &z).map_gens(|y| KhIGen::B(*y));
                    let qx = KhIChain::from(KhIGen::Q(*x));
                    let qtx = {
                        let tx = map(x);
                        KhIChain::from(KhIGen::Q(tx))
                    };
                    dx + qx + qtx
                },
                KhIGen::Q(x) => {
                    let z = KhChain::from(*x);
                    c.d(i, &z).map_gens(|y| KhIGen::Q(*y))
                }
            }
        };

        let inner = XChainComplex::new(summands, 1, move |i, z| { 
            z.apply(|x| d(i, x))
        });

        KhIComplex::new_impl(inner, vec![], deg_shift)
    }

    pub(crate) fn new_impl(inner: XChainComplex<KhIGen, R>, canon_cycles: Vec<KhIChain<R>>, deg_shift: (isize, isize)) -> Self { 
        Self { inner, canon_cycles, deg_shift }
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        let h0 = self.deg_shift.0;
        let h_min = self.support().min().unwrap_or(h0);
        let h_max = self.support().max().unwrap_or(h0);
        h_min ..= h_max
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        let q0 = self.deg_shift.1; 
        let q_itr = || self.support().flat_map(|i| self[i].gens().iter().map(|x| x.q_deg())); 
        let q_min = q_itr().min().unwrap_or(q0);
        let q_max = q_itr().max().unwrap_or(q0);
        q_min ..= q_max
    }

    pub fn canon_cycles(&self) -> &[KhIChain<R>] { 
        &self.canon_cycles
    }

    pub fn into_bigraded(self) -> XChainComplex2<KhIGen, R> {
        // TODO assert h == 0

        let h_range = self.h_range();
        let q_range = self.q_range().step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        let summands = Grid2::generate(support, |idx| { 
            let isize2(i, j) = idx;
            let gens = self[i].gens().iter().filter(|x| { 
                x.q_deg() == j
            }).cloned();
            XModStr::free(gens)
        });

        XChainComplex2::new(summands, isize2(1, 0), move |idx, x| { 
            let i = idx.0;
            let x = KhIChain::from(x.clone());
            let dx = self.d(i, &x);
            dx.into_iter().collect()
        })
    }

    pub fn homology(&self, with_trans: bool) -> KhIHomology<R>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
        let h = self.inner.reduced().homology(with_trans);
        KhIHomology::new_impl(h, self.h_range(), self.q_range())
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

#[cfg(all(test, feature = "old"))]
mod tests_old {
    #![allow(unused)]

    use itertools::Itertools;
    use yui::poly::HPoly;
    use yui::FF2;
    use num_traits::{Zero, One};
    use yui_homology::{ChainComplexCommon, ChainComplexTrait, DisplaySeq, DisplayTable, RModStr};
    use yui_link::Link;
    use super::*;

    #[test]
    fn complex_kh() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhIComplex::new_v1(&l, &h, false);

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
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::one();
        let c = KhIComplex::new_v1(&l, &h, false);

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
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        type P = HPoly<'H', R>;
        let h = P::variable();

        let c = KhIComplex::new_v1(&l, &h, false);

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
            [(1,5), (2,4)],
            Some(3)
        );

        type R = FF2;
        let h = R::zero();
        let c = KhIComplex::new_v1(&l, &h, true);

        assert_eq!(c.rank(0), 2);
        assert_eq!(c.rank(1), 5);
        assert_eq!(c.rank(2), 9);
        assert_eq!(c.rank(3), 10);
        assert_eq!(c.rank(4), 4);
        
        c.check_d_all();
    }

    #[test]
    fn complex_kh_bigr() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            None
        );

        type R = FF2;
        let h = R::zero();
        let c = KhIComplex::new_v1(&l, &h, false).into_bigraded();

        c.check_d_all();
    }

    #[test]
    fn complex_kh_red_bigr() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            Some(3)
        );

        type R = FF2;
        let h = R::zero();
        let c = KhIComplex::new_v1(&l, &h, true).into_bigraded();

        c.check_d_all();
    }

    #[test]
    fn canon_fbn() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            Some(3)
        );

        type R = FF2;
        let h = R::one();
        let c = KhIComplex::new_v1(&l, &h, false);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 4);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 0));
        assert!(zs[2].gens().all(|x| x.h_deg() == 1));
        assert!(zs[3].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = (i / 2) as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_fbn_red() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            Some(3)
        );

        type R = FF2;
        let h = R::one();
        let c = KhIComplex::new_v1(&l, &h, true);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 2);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = i as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_bn() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            Some(3)
        );

        type R = FF2;
        type P = HPoly<'H', R>;
        let h = P::variable();
        let c = KhIComplex::new_v1(&l, &h, false);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 4);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 0));
        assert!(zs[2].gens().all(|x| x.h_deg() == 1));
        assert!(zs[3].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = (i / 2) as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_bn_red() { 
        let l = InvLink::new(
            Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]), 
            [(1,5), (2,4)],
            Some(3)
        );

        type R = FF2;
        type P = HPoly<'H', R>;
        let h = P::variable();
        let c = KhIComplex::new_v1(&l, &h, true);
        
        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 2);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = i as isize;
            assert!(c.d(i, z).is_zero());
        }
    }
}