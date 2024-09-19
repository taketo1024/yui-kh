use std::collections::HashSet;
use std::marker::PhantomData;

use yui::{Ring, RingOps};
use yui_link::{Crossing, InvLink};

use crate::kh::v2::builder::TngComplexBuilder;
use crate::kh::KhComplex;

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    _r: PhantomData<R> 
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build(l: &InvLink, h: &R, t: &R, reduced: bool) { 
        assert!(l.link().is_knot(), "Only invertible knots are supported.");
        assert!(l.link().data().iter().all(|x| !x.is_resolved()));

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduced);
        let base_pt = if reduced { 
            // TODO check on-axis
            l.link().first_edge()
        } else { 
            None
        };

        let mut b_sym  = TngComplexBuilder::new(h, t, deg_shift, base_pt);
        let mut b_asym = TngComplexBuilder::new(h, t, (0, 0), None);

        let (x_sym, x_asym) = Self::split_crossings(l);
        b_sym.set_crossings(x_sym);
        b_sym.process();

        b_asym.set_crossings(x_asym);
        b_asym.process();

        let c_sym = b_sym.into_tng_complex();
        c_sym.print_d();

        let c_asym = b_asym.into_tng_complex();
        c_asym.print_d();

    }

    fn split_crossings(l: &InvLink) -> (HashSet<Crossing>, HashSet<Crossing>) { 
        let mut sym = HashSet::new();
        let mut asym = HashSet::new();
        let mut take = false; // a flag to collect only half of the asymmetric crossings.

        l.link().traverse_edges((0, 0), |i, j| { 
            let x = &l.link().data()[i];

            if i == l.inv_x(i) { // on-axis crossing
                // println!("{x} is sym.");

                if !sym.contains(x) {
                    // println!("  take {x}.");
                    sym.insert(x.clone());
                }
                
                take = !take;
                // println!("  passed on-axis crossing, take = {take}.");

            } else {  // off-axis crossing
                // println!("{x} is asym.");

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
    let l = InvLink::load("7_1").unwrap();
    SymTngBuilder::build(&l, &0, &0, false);
}