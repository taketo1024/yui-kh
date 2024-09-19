use std::collections::HashSet;
use std::marker::PhantomData;

use yui::{Ring, RingOps};
use yui_link::{Crossing, InvLink};

use crate::kh::v2::builder::TngComplexBuilder;
use crate::kh::v2::tng_complex::TngComplex;
use crate::kh::KhComplex;

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    _r: PhantomData<R> 
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_tng_complex(l: &InvLink, h: &R, t: &R, reduced: bool) -> TngComplex<R> { 
        assert!(l.link().is_knot(), "Only invertible knots are supported.");
        assert!(l.link().data().iter().all(|x| !x.is_resolved()));

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduced);
        let base_pt = if reduced { 
            // TODO check on-axis
            l.link().first_edge()
        } else { 
            None
        };

        let (x0, x1) = Self::split_crossings(l);

        let mut b0 = TngComplexBuilder::new(h, t, deg_shift, base_pt); // on-axis part
        b0.set_crossings(x0);
        b0.process();

        let c0 = b0.into_tng_complex();
        println!("center:\n{}", c0.desc_d());

        let mut b1 = TngComplexBuilder::new(h, t, (0, 0), None);       // half off-axis part
        b1.set_crossings(x1);
        b1.process();

        let c1 = b1.into_tng_complex();
        println!("left:\n{}", c1.desc_d());

        let c2 = c1.convert_edges(|e| l.inv_e(e));
        println!("right:\n{}", c2.desc_d());

        let c = c0.connect(&c1).connect(&c2);
        println!("combined:\n{}", c.desc_d());

        c
    }

    fn split_crossings(l: &InvLink) -> (HashSet<Crossing>, HashSet<Crossing>) { 
        let mut sym = HashSet::new();
        let mut asym = HashSet::new();
        let mut take = false; // a flag to collect only half of the asymmetric crossings.

        l.link().traverse_edges((0, 0), |i, j| { 
            let x = &l.link().data()[i];

            if i == l.inv_x(i) {
                // println!("{x} is on-axis.");

                take = !take;
                // println!("  passed on-axis crossing, take = {take}.");

                if !sym.contains(x) {
                    // println!("  take {x}.");
                    sym.insert(x.clone());
                }
            } else {
                // println!("{x} is off-axis.");

                let e = x.edge(j);
                if e == l.inv_e(e) { // on-axis edge
                    take = !take;
                    // println!("  passed on-axis edge, take = {take}.");
                }

                if take && !asym.contains(x) { 
                    // println!("  take {x}.");
                    asym.insert(x.clone());
                }
            }
        });

        assert_eq!(sym.len() + 2 * asym.len(), l.link().crossing_num());

        (sym, asym)
    }
}

#[test]
fn test() { 
    let l = InvLink::load("5_1").unwrap();
    SymTngBuilder::build_tng_complex(&l, &0, &0, false);
}