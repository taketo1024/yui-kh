use std::ops::RangeInclusive;

use yui::{EucRing, EucRingOps};
use yui_matrix::sparse::SpVec;

pub fn div_vec<R>(v: &SpVec<R>, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    v.iter().filter_map(|(_, a)|
        div(a, c)
    ).min()
}

fn div<R>(a: &R, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    if a.is_zero() { return None }

    let mut a = a.clone();
    let mut k = 0;

    while (&a % c).is_zero() { 
        a /= c;
        k += 1;
    }

    Some(k)
}

pub fn range_of<Idx, Itr>(itr: Itr) -> RangeInclusive<Idx>
where Idx: Ord + Default + Copy, Itr: IntoIterator<Item = Idx> { 
    let (min, max) = itr.into_iter().fold(None, |res, i| { 
        if let Some((min, max)) = res { 
            if i < min { 
                Some((i, max))
            } else if max < i { 
                Some((min, i))
            } else { 
                Some((min, max))
            }
        } else {
            Some((i, i))
        }
    }).unwrap_or((Idx::default(), Idx::default()));

    min ..= max
}