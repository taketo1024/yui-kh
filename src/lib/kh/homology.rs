use std::collections::HashMap;
use std::ops::{RangeInclusive, Index};
use delegate::delegate;
use cartesian::cartesian;

use yui_homology::{isize2, Grid2, GridTrait, RModStr, XHomology, XHomology2, XHomologySummand, XModStr};
use yui::{EucRing, EucRingOps};
use yui_link::Link;

use crate::kh::{KhGen, KhChainExt};

use super::{KhComplex, KhComplexBigraded};

pub struct KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: XHomology<KhGen, R>,
    h_range: RangeInclusive<isize>,
    q_range: RangeInclusive<isize>
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool, with_trans: bool) -> Self {
        Self::new_v2(l, h, t, reduced, with_trans)
    }

    pub fn new_v2(l: &Link, h: &R, t: &R, reduced: bool, with_trans: bool) -> Self {
        KhComplex::new_v2(l, h, t, reduced).homology(with_trans)
    }
    
    #[cfg(feature = "old")]
    pub fn new_v1(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        KhComplex::new_v1(l, h, t, reduced).homology(false)
    }

    pub(crate) fn new_impl(inner: XHomology<KhGen, R>, h_range: RangeInclusive<isize>, q_range: RangeInclusive<isize>) -> Self { 
        Self { inner, h_range, q_range }
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        self.h_range.clone()
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        self.q_range.clone()
    }

    pub fn inner(&self) -> &XHomology<KhGen, R> { 
        &self.inner
    }

    // TODO take self instead of &self
    pub fn into_bigraded(&self) -> KhHomologyBigraded<R> { 
        // TODO: check with_trans = true
        // TODO: check (h, t)

        let h_range = self.h_range();
        let q_range = self.q_range().step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        let table = self.collect_bigr_gens();

        let inner = Grid2::generate(support, move |idx| { 
            let i = idx.0;
            let Some(e) = table.get(&idx) else { 
                return XModStr::zero()
            };
            
            let (rank, tors, indices) = e;
            let gens = self[i].gens().clone(); 
            let trans = self[i].trans().map(|t|
                t.sub(indices)
            );
            XModStr::new(gens, *rank, tors.clone(), trans)
        });

        KhHomologyBigraded::new_impl(inner)
    }

    fn collect_bigr_gens(&self) -> HashMap<isize2, (usize, Vec<R>, Vec<usize>)> { 
        let mut table = HashMap::new();
        let init_entry = (0, vec![], vec![]);

        for i in self.support() { 
            let h = &self[i];
            let r = h.rank();
            let t = h.tors().len();

            for k in 0..r { 
                let z = h.gen_chain(k);
                let q = z.q_deg();
                let e = table.entry(isize2(i, q)).or_insert_with(|| init_entry.clone());
                e.0 += 1;
                e.2.push(k);
            }

            for k in 0..t { 
                let z = h.gen_chain(r + k);
                let q = z.q_deg();
                let e = table.entry(isize2(i, q)).or_insert_with(|| init_entry.clone());
                e.1.push(h.tors()[k].clone());
                e.2.push(r + k);
            }
        }

        table
    }
}

impl<R> GridTrait<isize> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Itr = std::vec::IntoIter<isize>;
    type Output = XHomologySummand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Output;
        }
    }
}

impl<R> Index<isize> for KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = XHomologySummand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

pub struct KhHomologyBigraded<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: XHomology2<KhGen, R>
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: Link, reduced: bool, with_trans: bool) -> Self {
        KhHomologyBigraded::new_v2(l, reduced, with_trans)
    }

    pub fn new_v2(l: Link, reduced: bool, with_trans: bool) -> Self {
        KhComplexBigraded::new_v2(l, reduced).homology(with_trans)
    }
    
    #[cfg(feature = "old")]
    pub fn new_v1(l: Link, reduced: bool, with_trans: bool) -> Self {
        KhComplexBigraded::new_v1(l, reduced).homology(with_trans)
    }

    pub(crate) fn new_impl(inner: XHomology2<KhGen, R>) -> Self { 
        Self { inner }
    }

    pub fn inner(&self) -> &XHomology2<KhGen, R> { 
        &self.inner
    }
}

impl<R> GridTrait<isize2> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Itr = std::vec::IntoIter<isize2>;
    type Output = XHomologySummand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize2) -> bool;
            fn get(&self, i: isize2) -> &Self::Output;
        }
    }
}

impl<R> Index<(isize, isize)> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = XHomologySummand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui::poly::HPoly;
    use yui::FF2;
    use yui_homology::RModStr;
    use yui_link::Link;
    use super::*;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::new_v2(&l, &0, &0, false, false);

        assert_eq!(h.h_range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::new_v2(&l, &0, &0, false, false);

        assert_eq!(h.h_range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::new_v2(&l, &0, &0, false, false);

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
        let h = KhHomology::new_v2(&l, &0, &0, false, false);

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
        let h = KhHomology::new_v2(&l, &0, &0, false, false);

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
        let h = KhHomologyBigraded::<i32>::new_v2(l, false, false);

        assert_eq!(h[(0,0)].rank(), 1);
        assert!(h[(0,0)].is_free());
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new_v2(l, false, false);

        assert_eq!(h[(0,-1)].rank(), 1);
        assert!(h[(0,-1)].is_free());
        assert_eq!(h[(0, 1)].rank(), 1);
        assert!(h[(0, 1)].is_free());
    }

    #[test]
    fn kh_unknot_bigr_red() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new_v2(l, true, false);

        assert_eq!(h[(0, 0)].rank(), 1);
        assert!(h[(0, 0)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new_v2(l, false, false);

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
        let h = KhHomologyBigraded::<i32>::new_v2(l, false, false);

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
        let h = KhHomologyBigraded::<i32>::new_v2(l, true, false);

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
        let h = KhHomologyBigraded::<i32>::new_v2(l, false, false);

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
        let h = KhHomologyBigraded::<i32>::new_v2(l, true, false);

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

   #[test]
   fn into_bigr() {
       let l = Link::trefoil();
       let (h, t) = (0, 0);
       let kh = KhHomology::new_v2(&l, &h, &t, false, true);
       let kh = kh.into_bigraded();

       assert_eq!(kh[(-3,-9)].rank(), 1);
       assert!(kh[(-3,-9)].is_free());
       assert_eq!(kh[(-2,-7)].rank(), 0);
       assert_eq!(kh[(-2,-7)].tors(), &vec![2]);
       assert_eq!(kh[(-2,-5)].rank(), 1);
       assert!(kh[(-2,-5)].is_free());
       assert_eq!(kh[( 0,-3)].rank(), 1);
       assert!(kh[( 0,-3)].is_free());
       assert_eq!(kh[( 0,-1)].rank(), 1);
       assert!(kh[( 0,-1)].is_free());
    }

    #[test]
    fn into_bigr_bn() {
        type R = FF2;
        type P = HPoly<'H', R>;

        let l = Link::trefoil();
        let (h, t) = (P::variable(), P::zero());
        let kh = KhHomology::new_v2(&l, &h, &t, false, true);
        let kh = kh.into_bigraded();

        assert_eq!(kh[(-2,-7)].rank(), 0);
        assert_eq!(kh[(-2,-7)].tors(), &vec![h.clone()]);
        assert_eq!(kh[(-2,-5)].rank(), 0);
        assert_eq!(kh[(-2,-5)].tors(), &vec![h.clone()]);
        assert_eq!(kh[( 0,-3)].rank(), 1);
        assert!(kh[( 0,-3)].is_free());
        assert_eq!(kh[( 0,-1)].rank(), 1);
        assert!(kh[( 0,-1)].is_free());
    }
 }

#[cfg(all(test, feature = "old"))]
mod tests_old {
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
        let h = KhHomologyBigraded::<i32>::new_v1(l, false, false);

        assert_eq!(h[(0,0)].rank(), 1);
        assert!(h[(0,0)].is_free());
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new_v1(l, false, false);

        assert_eq!(h[(0,-1)].rank(), 1);
        assert!(h[(0,-1)].is_free());
        assert_eq!(h[(0, 1)].rank(), 1);
        assert!(h[(0, 1)].is_free());
    }

    #[test]
    fn kh_unknot_bigr_red() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new_v1(l, true, false);

        assert_eq!(h[(0, 0)].rank(), 1);
        assert!(h[(0, 0)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new_v1(l, false, false);

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
        let h = KhHomologyBigraded::<i32>::new_v1(l, false, false);

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
        let h = KhHomologyBigraded::<i32>::new_v1(l, true, false);

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
        let h = KhHomologyBigraded::<i32>::new_v1(l, false, false);

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
        let h = KhHomologyBigraded::<i32>::new_v1(l, true, false);

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