#![allow(unused)]
use std::collections::{HashMap, HashSet};
use cartesian::cartesian;

use yui::bitseq::Bit;
use yui::{Ring, RingOps};
use yui_link::{Crossing, Edge, InvLink};

use crate::kh::v2::builder::TngComplexBuilder;
use crate::kh::v2::tng::{Tng, TngComp};
use crate::kh::v2::tng_complex::{TngComplex, TngKey};
use crate::kh::KhComplex;

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    e_map: HashMap<Edge, Edge>,
    key_map: HashMap<TngKey, TngKey>,
    on_axis: HashSet<Crossing>,
    off_axis: HashSet<Crossing>,
    complex: TngComplex<R>
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_tng_complex(l: &InvLink, h: &R, t: &R, reduced: bool, equivariant: bool) -> TngComplex<R> { 
        let mut b = Self::init(l, h, t, reduced);
        b.process_off_axis();
        if equivariant { 
            b.process_on_axis_equiv();
        } else { 
            b.process_on_axis_nonequiv();
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
        let key_map = HashMap::new();
        let (on_axis, off_axis) = Self::separate_crossings(l);

        SymTngBuilder { e_map, key_map, on_axis, off_axis, complex }
    }

    fn separate_crossings(l: &InvLink) -> (HashSet<Crossing>, HashSet<Crossing>) { 
        let mut on_axis = HashSet::new();
        let mut off_axis = HashSet::new();
        let mut take = false; // a flag to collect only half of the off-axis crossings.

        l.link().traverse_edges((0, 0), |i, j| { 
            let x = &l.link().data()[i];

            if i == l.inv_x(i) {
                take = !take;

                if !on_axis.contains(x) {
                    on_axis.insert(x.clone());
                }
            } else {
                let e = x.edge(j);
                if e == l.inv_e(e) { // on-axis edge
                    take = !take;
                }

                if take && !off_axis.contains(x) { 
                    off_axis.insert(x.clone());
                }
            }
        });

        assert_eq!(on_axis.len() + 2 * off_axis.len(), l.link().crossing_num());

        (on_axis, off_axis)
    }

    fn process_off_axis(&mut self) { 
        let (h, t) = self.complex.ht();

        let mut b = TngComplexBuilder::new(h, t, (0, 0), None); // half off-axis part
        b.set_crossings(self.off_axis.clone());
        b.process_all();

        let c1 = b.into_tng_complex();
        let c2 = c1.convert_edges(|e| self.inv_e(e));

        for ((k1, _), (k2, _)) in cartesian!(c1.iter_verts(), c2.iter_verts()) { 
            let k  = k1 + k2;
            let tk = k2 + k1;
            self.key_map.insert(k, tk);
        }

        self.complex.connect(c1);
        self.complex.connect(c2);
    }

    fn process_on_axis_nonequiv(&mut self) { 
        let complex = std::mem::take(&mut self.complex);
        
        let mut b = TngComplexBuilder::from(complex);
        b.set_crossings(self.on_axis.clone());
        b.process_all();

        self.complex = b.into_tng_complex();
    }

    fn process_on_axis_equiv(&mut self) { 
        let xs = self.on_axis.clone();

        for x in xs.iter() { 
            self.complex.append(x);
            self.update_key_map(x);

            // on-axis deloop
            // off-axis deloop
        }
    }

    fn update_key_map(&mut self, x: &Crossing) { 
        if x.is_resolved() { return }

        // (k <-> l) ==> (k0 <-> l0, k1 <-> l1) 
        self.key_map = self.key_map.iter().flat_map(|(k, l)| { 
            let mut e0 = (*k, *l);
            let mut e1 = e0;
            e0.0.state.push(Bit::Bit0);
            e0.1.state.push(Bit::Bit0);
            e1.0.state.push(Bit::Bit1);
            e1.1.state.push(Bit::Bit1);
            [e0, e1]
        }).collect();
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }

    fn inv_e(&self, e: Edge) -> Edge { 
        self.e_map[&e]
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
    
        // init_logger(log::LevelFilter::Trace);

        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (R::variable(), R::zero());
        let c = SymTngBuilder::build_tng_complex(&l, &h, &t, false, true);
        
        // c.print_d();
    }
}
