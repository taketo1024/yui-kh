#![allow(non_snake_case)]
use yui::poly::HPoly;
use yui::FF;
use yui_kh::{ssi_invariants, InvLink};

fn init_logger() { 
    use simplelog::*;
    TermLogger::init(
        LevelFilter::Debug,
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Always
    ).unwrap();
}

// cargo test -r -- --exact knotJ --nocapture --include-ignored
#[test]
#[ignore]
fn knotJ() { 
    type R = FF<2>;
    type P = HPoly<'H', R>;

    init_logger();
    
    let l = InvLink::knotJ();
    let c = P::variable();

    let ssi = ssi_invariants(&l, &c);

    assert_eq!(ssi.0, 0);
    assert_eq!(ssi.1, 2);
}

#[test]
#[ignore]
fn knotJm() { 
    type R = FF<2>;
    type P = HPoly<'H', R>;

    init_logger();
    
    let l = InvLink::knotJ().mirror();
    let c = P::variable();

    let ssi = ssi_invariants(&l, &c);
    
    assert_eq!(ssi.0, -2);
    assert_eq!(ssi.1, 0);
}