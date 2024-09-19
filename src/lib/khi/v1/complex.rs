use yui::{Ring, RingOps};
use yui_link::InvLink;

use crate::kh::KhComplex;
use crate::khi::v1::cube::KhICube;
use crate::khi::{KhIComplex, KhIGen};

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new_v1(l: &InvLink, h: &R, reduced: bool) -> Self { 
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        assert!(!reduced || l.base_pt().is_some());

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduced);
        let cube = KhICube::new(l, h, reduced, deg_shift);
        let inner = cube.into_complex();

        let canon_cycles = if l.base_pt().is_some() && l.link().is_knot() {
            let p = l.base_pt().unwrap();
            let zs = KhComplex::make_canon_cycles(l.link(), p, &R::zero(), h, reduced, deg_shift);
            Iterator::chain(
                zs.iter().map(|z| z.map_gens(|x| KhIGen::B(*x))),
                zs.iter().map(|z| z.map_gens(|x| KhIGen::Q(*x)))
            ).collect()
        } else { 
            vec![]
        };

        Self::new_impl(inner, canon_cycles, deg_shift)
    }
}