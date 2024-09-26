use std::collections::{HashMap, HashSet};

use yui::{Ring, RingOps};
use yui_link::{Crossing, Edge, InvLink};

use crate::kh::v2::builder::TngComplexBuilder;
use crate::kh::v2::tng::{Tng, TngComp};
use crate::kh::v2::tng_complex::{TngComplex, TngKey};
use crate::kh::KhComplex;

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    e_map: HashMap<Edge, Edge>,
    crossing_sym: HashSet<Crossing>,
    crossing_asym: HashSet<Crossing>,
    complex: TngComplex<R>
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_tng_complex(l: &InvLink, h: &R, t: &R, reduced: bool, equivariant: bool) -> TngComplex<R> { 
        let mut b = Self::init(l, h, t, reduced);
        b.process_asym();
        if equivariant { 
            b.process_sym_equiv();
        } else { 
            b.process_sym_nonequiv();
        }
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
        let e_map = l.link().edges().iter().map(|&e| (e, l.inv_e(e))).collect();
        let (sym, asym) = Self::split_crossings(l);

        SymTngBuilder { e_map, complex, crossing_sym: sym, crossing_asym: asym }
    }

    fn split_crossings(l: &InvLink) -> (HashSet<Crossing>, HashSet<Crossing>) { 
        let mut sym = HashSet::new();
        let mut asym = HashSet::new();
        let mut take = false; // a flag to collect only half of the asymmetric crossings.

        l.link().traverse_edges((0, 0), |i, j| { 
            let x = &l.link().data()[i];

            if i == l.inv_x(i) {
                take = !take;

                if !sym.contains(x) {
                    sym.insert(x.clone());
                }
            } else {
                let e = x.edge(j);
                if e == l.inv_e(e) { // on-axis edge
                    take = !take;
                }

                if take && !asym.contains(x) { 
                    asym.insert(x.clone());
                }
            }
        });

        assert_eq!(sym.len() + 2 * asym.len(), l.link().crossing_num());

        (sym, asym)
    }

    fn process_asym(&mut self) { 
        let (h, t) = self.complex.ht();

        let mut b = TngComplexBuilder::new(h, t, (0, 0), None); // half off-axis part
        b.set_crossings(self.crossing_asym.clone());
        b.process_all();

        let c1 = b.into_tng_complex();
        let c2 = c1.convert_edges(|e| self.inv_e(e));

        self.complex.connect(c1);
        self.complex.connect(c2);
    }

    fn process_sym_nonequiv(&mut self) { 
        let complex = std::mem::take(&mut self.complex);
        
        let mut b = TngComplexBuilder::from(complex);
        b.set_crossings(self.crossing_sym.clone());
        b.process_all();

        self.complex = b.into_tng_complex();
    }

    fn process_sym_equiv(&mut self) { 
        for x in self.crossing_sym.iter() { 
            self.complex.append(x);

            // on-axis symmetric deloop
            while let Some((k, r)) = self.complex.find_comp(|c| 
                self.is_symmetric_circle(c)
            ) { 
                let ks = self.complex.deloop(&k, r);
                for k in ks { 
                    if let Some((i, j)) = self.complex.find_inv_edge_with(&k, |i, j| {
                        self.is_symmetric_cob(i, j)
                    }) { 
                        self.complex.eliminate(&i, &j);
                    }
                }
            }

            // off-axis equivariant deloop
        }
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }

    fn inv_e(&self, e: Edge) -> Edge { 
        self.e_map[&e]
    }

    fn is_symmetric_circle(&self, c: &TngComp) -> bool { 
        c.is_circle() && &c.convert_edges(|e| self.inv_e(e)) == c
    }

    fn is_symmetric_tng(&self, c: &Tng) -> bool { 
        &c.convert_edges(|e| self.inv_e(e)) == c
    }

    fn is_symmetric_cob(&self, i: &TngKey, j: &TngKey) -> bool { 
        let t1 = self.complex.vertex(i).tng();
        let t2 = self.complex.vertex(j).tng();
        self.is_symmetric_tng(t1) && self.is_symmetric_tng(t2) 
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests { 
    use super::*;
    use num_traits::Zero;
    use yui_homology::{ChainComplexCommon, DisplaySeq};

    fn init_logger(l: log::LevelFilter) {
        use simplelog::*;
        TermLogger::init(
            l,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto
        ).unwrap()
    }

    #[test]
    fn test() { 
        use yui::FF2;
        use yui::poly::HPoly;
        type R = HPoly<'H', FF2>;
    
        init_logger(log::LevelFilter::Trace);

        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (R::variable(), R::zero());
        let c = SymTngBuilder::build_tng_complex(&l, &h, &t, false, true);
        
        c.print_d();
    }
}
