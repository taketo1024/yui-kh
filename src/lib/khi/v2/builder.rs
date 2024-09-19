use std::marker::PhantomData;

use yui::{Ring, RingOps};
use yui_link::InvLink;

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    _r: PhantomData<R> 
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build(l: &InvLink, h: &R, t: &R, reduced: bool) { 
        todo!()
    }
}

#[test]
fn test() { 
    let l = InvLink::load("5_1").unwrap();
    SymTngBuilder::build(&l, &0, &0, false);
}