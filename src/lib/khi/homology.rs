use std::collections::HashMap;
use std::ops::{Index, RangeInclusive};
use cartesian::cartesian;
use delegate::delegate;
use yui::{EucRing, EucRingOps};
use yui_homology::{isize2, Grid2, GridTrait, RModStr, SimpleRModStr, XHomology, XHomologySummand};
use yui_link::InvLink;
use crate::khi::{KhIChainExt, KhIComplex, KhIGen};

pub struct KhIHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: XHomology<KhIGen, R>,
    h_range: RangeInclusive<isize>,
    q_range: RangeInclusive<isize>
}

impl<R> KhIHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &InvLink, h: &R, reduced: bool, with_trans: bool) -> Self {
        KhIComplex::new(l, h, reduced).homology(with_trans)
    }

    pub fn new_impl(inner: XHomology<KhIGen, R>, h_range: RangeInclusive<isize>, q_range: RangeInclusive<isize>) -> Self {
        Self { inner, h_range, q_range }
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        self.h_range.clone()
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        self.q_range.clone()
    }

    pub fn inner(&self) -> &XHomology<KhIGen, R> { 
        &self.inner
    }

    pub fn gen_table(&self) -> Grid2<SimpleRModStr<R>> { 
        // TODO: check with_trans = true

        let h_range = self.h_range();
        let q_range = self.q_range().step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        let mut table: HashMap<isize2, (usize, Vec<R>)> = HashMap::new();
        let e = (0, vec![]);

        for i in self.support() { 
            let h = &self[i];
            let r = h.rank();
            let t = h.tors().len();

            for k in 0..r { 
                let z = h.gen_chain(k);
                let q = z.q_deg();
                let e = table.entry(isize2(i, q)).or_insert_with(|| e.clone());
                e.0 += 1;
            }

            for k in 0..t { 
                let z = h.gen_chain(r + k);
                let q = z.q_deg();
                let e = table.entry(isize2(i, q)).or_insert_with(|| e.clone());
                e.1.push(h.tors()[k].clone());
            }
        }

        Grid2::generate(support, move |idx| { 
            let Some(e) = table.remove(&idx) else { 
                return SimpleRModStr::zero()
            };
            SimpleRModStr::new(e.0, e.1, None)
        })
    }
}

impl<R> GridTrait<isize> for KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Itr = std::vec::IntoIter<isize>;
    type Output = XHomologySummand<KhIGen, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Output;
        }
    }
}

impl<R> Index<isize> for KhIHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = XHomologySummand<KhIGen, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}
