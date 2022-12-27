use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Display;
use std::ops::{Add, Sub, Index, IndexMut};
use std::vec::IntoIter;
use crate::math::traits::{MathElem, Ring, RingOps};
use crate::utils::format::superscript;

pub trait RModStr: Sized + Display
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> { 
    type R;

    fn zero() -> Self;
    fn rank(&self) -> usize;
    fn tors(&self) -> &Vec<Self::R>;

    fn is_free(&self) -> bool { 
        self.tors().is_empty()
    }

    fn is_zero(&self) -> bool { 
        self.rank() == 0 && self.is_free()
    }

    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use grouping_by::GroupingBy;
        let rank = self.rank();
        let tors = self.tors();

        if rank == 0 && tors.is_empty() { 
            return f.write_str("0")
        }
    
        let mut res = vec![];
        let symbol = Self::R::symbol();
    
        if rank > 1 {
            let str = format!("{}{}", symbol, superscript(rank as isize));
            res.push(str);
        } else if rank == 1 { 
            let str = format!("{}", symbol);
            res.push(str);
        }
    
        for (t, r) in tors.iter().counter(|&t| t) { 
            let str = if r > 1 { 
                format!("({}/{}){}", symbol, t, superscript(r as isize))
            } else { 
                format!("({}/{})", symbol, t)
            };
            res.push(str);
        }
    
        let str = res.join(" ⊕ ");
        f.write_str(&str) 
    }
}

#[derive(Debug, Clone)]
pub struct GenericRModStr<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize,
    tors: Vec<R>
}

impl<R> GenericRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(rank: usize, tors: Vec<R>) -> Self { 
        Self { rank, tors }
    }
}

impl<R> RModStr for GenericRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn zero() -> Self { 
        Self { rank: 0, tors: vec![] }
    }

    fn rank(&self) -> usize { 
        self.rank
    }

    fn tors(&self) -> &Vec<R> {
        &self.tors
    }
}

impl<R> Display for GenericRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
}

pub trait AdditiveIndex: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>
{}

impl <T> AdditiveIndex for T
where T: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>{}

pub trait AdditiveIndexRange: Iterator + Clone 
where Self::Item: AdditiveIndex
{}

impl <T> AdditiveIndexRange for T
where T: Iterator + Clone, T::Item: AdditiveIndex
{}

pub trait GradedRModStr: Index<Self::Index>
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
    Self::Index: AdditiveIndex,
    Self::IndexRange: AdditiveIndexRange<Item = Self::Index>
{
    type R;
    type Index;
    type IndexRange;

    fn in_range(&self, k: Self::Index) -> bool;
    fn range(&self) -> Self::IndexRange;

    fn is_free(&self) -> bool { 
        self.range().all(|i| self[i].is_free())
    }

    fn is_zero(&self) -> bool { 
        self.range().all(|i| self[i].is_zero())
    }
}

pub struct RModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    range: I,
    grid: HashMap<I::Item, S>,
    zero: S
}

impl<S, I> RModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    pub fn new<F>(range: I, mut f: F) -> Self
    where F: FnMut(I::Item) -> S {
        let grid = range.clone().filter_map(|i| {
            let s = f(i);
            if !s.is_zero() { Some((i, s)) } else { None }
        }).collect();
        let zero = S::zero();
        Self { range, grid, zero }
    }
}

impl<S, I> Index<I::Item> for RModGrid<S, I>
where
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type Output = S;

    fn index(&self, index: I::Item) -> &Self::Output {
        if let Some(s) = self.grid.get(&index) { 
            s
        } else { 
            &self.zero
        }
    }
}

impl<S, I> IndexMut<I::Item> for RModGrid<S, I>
where
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    fn index_mut(&mut self, index: I::Item) -> &mut Self::Output {
        self.grid.get_mut(&index).unwrap()
    }
}

impl<S, I> GradedRModStr for RModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type R = S::R;
    type Index = I::Item;
    type IndexRange = I;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.contains_key(&k)
    }

    fn range(&self) -> Self::IndexRange {
        self.range.clone()
    }
}


use derive_more::{Display, Add, Sub};
use itertools::Itertools;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Display, Add, Sub)]
#[display(fmt = "({}, {})", _0, _1)]
pub struct Idx2(pub isize, pub isize);

impl Idx2 { 
    pub fn iterate(from: Idx2, to: Idx2, step:(usize, usize)) -> Idx2Range {
        (from.1 ..= to.1).step_by(step.1).flat_map(|j| { 
            (from.0 ..= to.0).step_by(step.0).map(move |i| Idx2(i, j))
        }).collect_vec().into_iter()
    }

    pub fn as_tuple(&self) -> (isize, isize) {
        (self.0, self.1)
    }
}

impl From<[isize; 2]> for Idx2 { 
    fn from(idx: [isize; 2]) -> Self {
        Self(idx[0], idx[1])
    }
}

pub type Idx2Range = IntoIter<Idx2>;