use std::collections::HashSet;

use yui::{Ring, RingOps};
use yui_link::{Crossing, InvLink};

use crate::kh::v2::builder::TngComplexBuilder;
use crate::kh::v2::tng_complex::TngComplex;
use crate::kh::KhComplex;

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    link: InvLink,
    crossing_sym: HashSet<Crossing>,
    crossing_asym: HashSet<Crossing>,
    complex: TngComplex<R>
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_tng_complex(l: &InvLink, h: &R, t: &R, reduced: bool) -> TngComplex<R> { 
        let mut b = Self::init(l, h, t, reduced);
        b.process_asym();
        b.process_sym();
        b.into_tng_complex()
    }

    pub fn init(l: &InvLink, h: &R, t: &R, reduced: bool) -> SymTngBuilder<R> { 
        assert!(l.link().is_knot(), "Only invertible knots are supported.");
        assert!(l.link().data().iter().all(|x| !x.is_resolved()));

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduced);
        let base_pt = if reduced { 
            // TODO check on-axis
            l.link().first_edge()
        } else { 
            None
        };

        let complex = TngComplex::init(h, t, deg_shift, base_pt);
        let (sym, asym) = Self::split_crossings(l);

        SymTngBuilder { link: l.clone(), complex, crossing_sym: sym, crossing_asym: asym }
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

    fn process_asym(&mut self) { 
        let (h, t) = self.complex.ht();
        let l = &self.link;
        let x_asym = std::mem::take(&mut self.crossing_asym);

        let mut b = TngComplexBuilder::new(h, t, (0, 0), None); // half off-axis part
        b.set_crossings(x_asym);
        b.process_all();

        let c1 = b.into_tng_complex();
        let c2 = c1.convert_edges(|e| l.inv_e(e));

        let c = self.complex.connect(&c1).connect(&c2);
        self.complex = c;
    }

    fn process_sym(&mut self) { 
        let complex = std::mem::take(&mut self.complex);
        let x_sym = std::mem::take(&mut self.crossing_sym);

        let mut b = TngComplexBuilder::from(complex);
        b.set_crossings(x_sym);
        b.process_all();

        self.complex = b.into_tng_complex();
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests { 
    use super::*;
    use num_traits::Zero;
    use yui_homology::{ChainComplexCommon, DisplaySeq};

    #[test]
    fn test() { 
        use yui::FF2;
        use yui::poly::HPoly;
        type R = HPoly<'H', FF2>;
    
        let l = InvLink::sinv_knot_from_code([
            [1,27,2,26],[19,2,20,3],[3,13,4,12],[4,31,5,32],[30,5,31,6],
            [13,7,14,6],[8,27,9,28],[9,1,10,34],[10,18,11,17],[24,11,25,12],
            [14,21,15,22],[28,16,29,15],[33,16,34,17],[18,26,19,25],[20,8,21,7],
            [29,23,30,22],[23,33,24,32]
        ]);
        let (h, t) = (R::variable(), R::zero());
        let c = SymTngBuilder::build_tng_complex(&l, &h, &t, false);
        
        c.print_d();
    }
}
