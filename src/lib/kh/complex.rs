use std::collections::HashMap;
use std::ops::{RangeInclusive, Index};
use cartesian::cartesian;

use delegate::delegate;
use yui::lc::Lc;
use yui::{Ring, RingOps, EucRing, EucRingOps};
use yui_link::Link;
use yui_homology::{isize2, ChainComplexTrait, Grid2, GridTrait, RModStr, SimpleRModStr, XChainComplex, XChainComplex2, XChainComplexSummand, XModStr};
use yui_matrix::sparse::SpMat;

use crate::kh::{KhGen, KhHomology, KhHomologyBigraded};

pub type KhChain<R> = Lc<KhGen, R>;
pub trait KhChainExt { 
    fn h_deg(&self) -> isize;
    fn q_deg(&self) -> isize;
}

impl<R> KhChainExt for KhChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn h_deg(&self) -> isize {
        self.gens().map(|x| x.h_deg()).min().unwrap_or(0)
    }
    
    fn q_deg(&self) -> isize {
        self.gens().map(|x| x.q_deg()).min().unwrap_or(0)
    }
}

pub type KhComplexSummand<R> = XChainComplexSummand<KhGen, R>;

pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: XChainComplex<KhGen, R>,
    canon_cycles: Vec<KhChain<R>>,
    reduced: bool,
    deg_shift: (isize, isize)
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(link: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        Self::new_v2(link, h, t, reduced)
    }

    pub(crate) fn new_impl(inner: XChainComplex<KhGen, R>, canon_cycles: Vec<KhChain<R>>, reduced: bool, deg_shift: (isize, isize)) -> Self { 
        KhComplex { inner, canon_cycles, reduced, deg_shift }
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
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

    pub fn inner(&self) -> &XChainComplex<KhGen, R> {
        &self.inner
    }

    pub fn gen_table(&self) -> Grid2<SimpleRModStr<R>> { 
        let h_range = self.h_range();
        let q_range = self.q_range().step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        let mut table: HashMap<isize2, (usize, Vec<R>)> = HashMap::new();
        let e = (0, vec![]);

        for i in self.support() { 
            let c = &self[i];
            let r = c.rank();

            for k in 0..r { 
                let z = c.gens()[k];
                let q = z.q_deg();
                let e = table.entry(isize2(i, q)).or_insert_with(|| e.clone());
                e.0 += 1;
            }
        }

        Grid2::generate(support, move |idx| { 
            let Some(e) = table.remove(&idx) else { 
                return SimpleRModStr::zero()
            };
            SimpleRModStr::new(e.0, e.1, None)
        })
    }

    pub fn into_bigraded(self) -> KhComplexBigraded<R> {
        let reduced = self.reduced;

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

        let canon_cycles = self.canon_cycles.clone();

        let inner = XChainComplex2::new(summands, isize2(1, 0), move |idx, x| { 
            let i = idx.0;
            let x = KhChain::from(x.clone());
            let dx = self.d(i, &x);
            dx.into_iter().collect()
        });

        KhComplexBigraded { inner, canon_cycles, reduced }
    }

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.signed_crossing_nums();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg;
        let e = if reduced { 1 } else { 0 };
        (h, q + e)
    }

    delegate! { 
        to self.inner { 
            pub fn d(&self, i: isize, z: &KhChain<R>) -> KhChain<R>;
        }
    }
}

impl<R> Index<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner {
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<isize>;
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Output;
        }
    }
}

impl<R> ChainComplexTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Element = KhChain<R>;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize) -> usize;
            fn d_deg(&self) -> isize;
            fn d(&self, i: isize, z: &Self::Element) -> Self::Element;
            fn d_matrix(&self, i: isize) -> SpMat<R>;
        }
    }
}

impl<R> KhComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self, with_trans: bool) -> KhHomology<R> {
        let h = self.inner.homology(with_trans);
        KhHomology::new_impl(h, self.h_range(), self.q_range())
    }
}

pub struct KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: XChainComplex2<KhGen, R>,
    canon_cycles: Vec<KhChain<R>>,
    reduced: bool,
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: Link, reduced: bool) -> Self { 
        Self::new_v2(l, reduced)
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn inner(&self) -> &XChainComplex2<KhGen, R> {
        &self.inner
    }
}

impl<R> Index<(isize, isize)> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<isize2>;
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize2) -> bool;
            fn get(&self, i: isize2) -> &Self::Output;
        }
    }
}

impl<R> ChainComplexTrait<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Element = KhChain<R>;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize2) -> usize;
            fn d_deg(&self) -> isize2;
            fn d(&self, i: isize2, z: &Self::Element) -> Self::Element;
            fn d_matrix(&self, i: isize2) -> SpMat<Self::R>;
        }
    }
}

impl<R> KhComplexBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self, with_trans: bool) -> KhHomologyBigraded<R> {
        let h = self.inner.homology(with_trans);
        KhHomologyBigraded::new_impl(h)
    }
}