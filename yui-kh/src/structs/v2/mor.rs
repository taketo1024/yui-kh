use yui::{Ring, RingOps};
use yui::lc::LinComb;

use super::cob::{Cob, Dot, Bottom, CobComp};
use super::tng::{Tng, TngComp};

pub type Mor<R> = LinComb<Cob, R>; // R-linear combination of cobordisms.

pub trait MorTrait: Sized {
    type R;
    fn src(&self) -> Tng;
    fn tgt(&self) -> Tng;
    fn is_closed(&self) -> bool;
    fn is_invertible(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn map_cob<F>(self, f: F) -> Self where F: Fn(&mut Cob);
    fn connect(self, c: &Cob) -> Self;
    fn connect_comp(self, c: &CobComp) -> Self;
    fn cap_off(self, b: Bottom, c: &TngComp, dot: Dot) -> Self;
    fn part_eval(self, h: &Self::R, t: &Self::R) -> Self;
    fn eval(&self, h: &Self::R, t: &Self::R) -> Self::R;
}

impl<R> MorTrait for Mor<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn src(&self) -> Tng {
        let Some((c, _)) = self.iter().next() else { 
            return Tng::empty()
        };
        c.src()
    }

    fn tgt(&self) -> Tng { 
        let Some((c, _)) = self.iter().next() else { 
            return Tng::empty()
        };
        c.tgt()
    }

    fn is_closed(&self) -> bool { 
        self.iter().all(|(f, _)| f.is_closed())
    }

    fn is_invertible(&self) -> bool { 
        self.len() == 1 && 
        self.iter().next().map(|(c, a)| 
            c.is_invertible() && a.is_unit()
        ).unwrap_or(false)
    }

    fn inv(&self) -> Option<Self> { 
        if let Some((Some(cinv), Some(ainv))) = self.iter().next().map(|(c, a)| 
            (c.inv(), a.inv())
        ) { 
            let inv = Mor::from((cinv, ainv));
            Some(inv)
        } else { 
            None
        }
    }

    fn map_cob<F>(self, f: F) -> Self 
    where F: Fn(&mut Cob) {
        self.into_map(|mut cob, r| { 
            f(&mut cob);
            if cob.is_zero() { 
                (cob, R::zero())
            } else { 
                (cob, r)
            }
        })
    }

    fn connect(self, c: &Cob) -> Self {
        self.map_cob(|cob| cob.connect(c.clone()) )
    }

    fn connect_comp(self, c: &CobComp) -> Self {
        self.map_cob(|cob| cob.connect_comp(c.clone()) )
    }

    fn cap_off(self, b: Bottom, c: &TngComp, dot: Dot) -> Self {
        self.map_cob(|cob| cob.cap_off(b, c, dot) )
    }

    fn part_eval(self, h: &Self::R, t: &Self::R) -> Self {
        if self.keys().any(|c| c.should_part_eval()) { 
            self.into_iter().map(|(cob, r)|
                cob.part_eval(h, t) * r
            ).sum()
        } else { 
            self
        }
    }

    fn eval(&self, h: &R, t: &R) -> R {
        self.iter().map(|(c, a)| { 
            a * c.eval(h, t)
        }).sum()
    }
}