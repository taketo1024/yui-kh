use std::collections::HashMap;
use std::mem;

use itertools::{Itertools, Either};
use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3};
use yui_matrix::sparse::{SpMat, SpVec, MatType, Trans};

use crate::ReducedComplexBase;

use super::graded::Graded;
use super::utils::{ChainReducer, HomologyCalc};
use super::homology::{HomologySummand, HomologyBase};

pub type ChainComplex<R>  = ChainComplexBase<isize,  R>;
pub type ChainComplex2<R> = ChainComplexBase<isize2, R>;
pub type ChainComplex3<R> = ChainComplexBase<isize3, R>;

pub trait ChainComplexTrait<I>: Graded<I> + Sized
where I: Deg, Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> { 
    type R;

    fn d_deg(&self) -> I;
    fn d_matrix(&self, i: I) -> &SpMat<Self::R>;

    fn rank(&self, i: I) -> usize { 
        if self.is_supported(i) { 
            self.d_matrix(i).shape().1 // column
        } else { 
            0
        }
    }

    fn d(&self, i: I, v: &SpVec<Self::R>) -> SpVec<Self::R> {
        assert_eq!(self.rank(i), v.dim());
        let d = self.d_matrix(i);
        d * v
    }

    fn is_cycle(&self, i: I, v: &SpVec<Self::R>) -> bool { 
        self.d(i, v).is_zero()
    }

    fn check_d_at(&self, i: I) { 
        let i1 = i + self.d_deg();
        if !(self.is_supported(i) && self.is_supported(i1)) {
            return 
        }

        let d0 = self.d_matrix(i);
        let d1 = self.d_matrix(i1);
        let res = d1 * d0;

        assert!( res.is_zero(), "d² is non-zero at {i}." );
    }

    fn check_d_all(&self) {
        for i in self.support() { 
            self.check_d_at(i);
        }
    }
    
    fn display_d_at(&self, i: I) -> String {
        let c = |i| self.display_at(i).unwrap();
        let c0 = c(i);
        let c1 = c(i + self.d_deg());
        let d = self.d_matrix(i).to_dense();

        if d.is_zero() { 
            format!("C[{i}] {c0} -> {c1}; zero.")
        } else { 
            format!("C[{i}] {c0} -> {c1}\n{d}.")
        }
    }

    fn display_d(&self) -> String { 
        self.support().filter_map(|i| 
            if self.rank(i) > 0 && self.rank(i + self.d_deg()) > 0 {
                Some(self.display_d_at(i))
            } else { 
                None
            }
        ).join("\n\n")
    }

    fn print_d(&self) {
        println!("{}", self.display_d());
    }

    fn reduced(&self, with_trans: bool) -> ReducedComplexBase<I, Self::R> { 
        ChainReducer::reduce(self, with_trans)
    }

    fn reduced_by<F>(&self, trans: F) -> ReducedComplexBase<I, Self::R>
    where F: FnMut(I) -> Trans<Self::R> {
        ReducedComplexBase::reduced_by(self, trans)
    }
}

pub struct ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    support: Vec<I>,
    d_deg: I,
    d_matrix: HashMap<I, SpMat<R>>,
    zero_d: SpMat<R>
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F>(support: It, d_deg: I, mut d_matrix: F) -> Self
    where 
        It: Iterator<Item = I>, 
        F: FnMut(I) -> SpMat<R>
    {
        let support = support.collect_vec();
        let mut d_matrix = HashMap::from_iter( 
            support.iter().map(|&i| (i, d_matrix(i))) 
        );

        // add zero maps
        for &i in support.iter() {
            if !d_matrix.contains_key(&(i - d_deg)) { 
                let r = d_matrix[&i].shape().1;
                let d = SpMat::zero((r, 0));
                d_matrix.insert(i - d_deg, d);
            }
        }
        let zero_d = SpMat::zero((0, 0));

        Self { 
            support, d_deg, d_matrix, zero_d
        }
    }
}

impl<I, R> Graded<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<I>;

    fn support(&self) -> Self::Itr {
        self.support.clone().into_iter()
    }

    fn display_at(&self, i: I) -> Option<String> {
        use yui_utils::superscript;

        let symbol = R::math_symbol();
        let rank = self.rank(i);
        if rank > 1 {
            let f = format!("{}{}", symbol, superscript(rank as isize));
            Some(f)
        } else if rank == 1 { 
            Some(symbol)
        } else { 
            None
        }
    }
}

impl<I, R> ChainComplexTrait<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn d_deg(&self) -> I { 
        self.d_deg
    }

    fn d_matrix(&self, i: I) -> &SpMat<R> { 
        if let Some(d) = self.d_matrix.get(&i) {
            d
        } else { 
            &self.zero_d
        }
    }
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology_at(&self, i: I, with_trans: bool) -> HomologySummand<R> {
        let i0 = i - self.d_deg();
        let d0 = self.d_matrix(i0);
        let d1 = self.d_matrix(i);
        HomologyCalc::calculate(d0, d1, with_trans)
    }

    pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R> {
        HomologyBase::new(&self, with_trans)
    }
}

impl<R> ChainComplexBase<isize, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_mats(d_deg: isize, offset: isize, mut mats: Vec<SpMat<R>>) -> Self { 
        let n = mats.len() as isize;
        let range = offset .. offset + n;
        let range = if d_deg.is_positive() { 
            Either::Left(range) 
        } else { 
            Either::Right(range.rev())
        };
        Self::new(
            range, d_deg, 
            move |i| {
                let i = (i - offset) as usize;
                mem::take(&mut mats[i])
            }
        )
    }
}

#[cfg(test)]
pub(crate) mod tests { 
    use super::*;

    #[test]
    fn d3() { 
        let c = ChainComplex::<i64>::d3();

        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);
        assert_eq!(c.rank(3), 1);

        c.check_d_all();
    }

    #[test]
    fn s2() { 
        let c = ChainComplex::<i64>::s2();

        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);
        assert_eq!(c.rank(3), 0);

        c.check_d_all();
    }

    #[test]
    fn t2() { 
        let c = ChainComplex::<i64>::t2();

        assert_eq!(c.rank(0), 9);
        assert_eq!(c.rank(1), 27);
        assert_eq!(c.rank(2), 18);
        assert_eq!(c.rank(3), 0);

        c.check_d_all();
    }

    #[test]
    fn rp2() { 
        let c = ChainComplex::<i64>::rp2();
        
        assert_eq!(c.rank(0), 6);
        assert_eq!(c.rank(1), 15);
        assert_eq!(c.rank(2), 10);
        assert_eq!(c.rank(3), 0);

        c.check_d_all();
    }
}