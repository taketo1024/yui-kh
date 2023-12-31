use std::collections::HashMap;
use std::fmt::Display;
use itertools::Itertools;
use num_traits::Zero;
use yui::{Ring, RingOps};
use yui_link::{State, Crossing};
use yui::macros::map;

use crate::{KhEnhState, KhAlgGen, KhChain};

use super::cob::{Cob, Bottom, Dot};
use super::mor::{Mor, MorTrait};
use super::tng::{TngComp, Tng};

// element in C as a cobordism ∅ → C.
pub struct TngElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    state: State,
    value: Cob,                        // precomposed at the final step.
    mors: HashMap<KhEnhState, Mor<R>>, // src must partially match init_cob. 
    x_count: usize
}

impl<R> TngElem<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn init(state: State, value: Cob) -> Self { 
        let f = Mor::from(Cob::empty());
        let mors = map! { KhEnhState::init() => f };
        let x_count = 0;

        Self{ state, value, mors, x_count }
    }

    pub fn append(&mut self, x: &Crossing) { 
        if !x.is_resolved() { 
            self.append_x(x)
        } else { 
            self.append_a(x)
        }
    }

    fn append_x(&mut self, x: &Crossing) {
        assert!(!x.is_resolved());

        let i = self.x_count;
        let r = self.state[i];

        let a = x.resolved(r);
        let tng = Tng::from_a(&a);
        let id = Cob::id(&tng);

        let mors = std::mem::take(&mut self.mors);
        self.mors = mors.into_iter().map(|(mut k, f)| {
            k.state.push(r);
            let f = f.connect(&id);
            (k, f)
        }).collect();

        self.x_count += 1;
    }

    fn append_a(&mut self, x: &Crossing) {
        assert!(x.is_resolved());

        let tng = Tng::from_a(x);
        let id = Cob::id(&tng);

        let mors = std::mem::take(&mut self.mors);
        self.mors = mors.into_iter().map(|(k, f)| {
            let f = f.connect(&id);
            (k, f)
        }).collect();
    }

    pub fn deloop(&mut self, k: &KhEnhState, c: &TngComp, reduced: bool) {
        let Some(f) = self.mors.remove(k) else { return };

        let (k0, f0) = self.deloop_for(k, &f, c, KhAlgGen::X, Dot::None);
        self.mors.insert(k0, f0);

        if !reduced { 
            let (k1, f1) = self.deloop_for(k, &f, c, KhAlgGen::I, Dot::Y);
            self.mors.insert(k1, f1);    
        }
    }

    fn deloop_for(&self, k: &KhEnhState, f: &Mor<R>, c: &TngComp, label: KhAlgGen, dot: Dot) -> (KhEnhState, Mor<R>) { 
        let mut k_new = *k;
        k_new.label.push(label);

        let f_new = f.clone().cap_off(Bottom::Tgt, c, dot);
        (k_new, f_new)
    }

    pub fn eliminate(&mut self, i: &KhEnhState, j: &KhEnhState, i_out: &HashMap<KhEnhState, Mor<R>>) {
        // mors into i can be simply dropped.
        self.mors.remove(i);

        // mors into j must be redirected by -ca^{-1}
        let Some(f) = self.mors.remove(j) else { return };

        let a = &i_out[&j];
        let ainv = a.inv().unwrap();

        for (k, c) in i_out.iter() { 
            if k == j { continue }

            let caf = c * &ainv * &f;
            let d = if let Some(d) = self.mors.get(k) { 
                d - caf
            } else { 
                -caf
            };

            if !d.is_zero() { 
                self.mors.insert(*k, d);
            } else { 
                self.mors.remove(k);
            }
        }
    }

    pub fn finalize(&mut self) { 
        let val = std::mem::take(&mut self.value);
        let mut mors = std::mem::take(&mut self.mors);

        mors = mors.into_iter().map(|(k, f)|
            (k, f.map_cob(|c| *c = &*c * &val))
        ).collect();
        mors.retain(|_, f| !f.is_zero());

        self.mors = mors;
    }

    pub fn eval(&self, h: &R, t: &R) -> KhChain<R> {
        assert!(self.value.is_empty());
        assert!(self.mors.values().all(|f| f.is_closed()));

        KhChain::from_iter(self.mors.iter().map(|(k, f)|
            (*k, f.eval(h, t))
        ))
    }
}

impl<R> Display for TngElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mors = self.mors.iter().map(|(k, f)| { 
            format!("{}: {}", k, f)
        }).join(", ");
        write!(f, "[{}]", mors)
    }
}