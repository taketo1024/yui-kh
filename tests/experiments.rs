#![allow(non_snake_case)]
use itertools::Itertools;
use yui::poly::HPoly;
use yui::FF2;
use yui_link::InvLink;
use yui_kh::ssi_invariants;

type R = FF2;
type P = HPoly<'H', R>;

fn init_logger() { 
    use simplelog::*;
    TermLogger::init(
        LevelFilter::Debug,
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Always
    ).unwrap();
}

#[test]
fn shit() { 
    let l = [
        [0,26,1,25],[18,1,19,2],[2,12,3,11],[3,30,4,31],[29,4,30,5],
        [12,6,13,5],[7,26,8,27],[8,0,9,33],[9,17,10,16],[23,10,24,11],
        [13,20,14,21],[27,15,28,14],[32,15,33,16],[17,25,18,24],[19,7,20,6],
        [28,22,29,21],[22,32,23,31]
    ];
    let l = format!("[{}]", l.iter().map(|x| format!("[{}]", x.map(|i| (i + 1).to_string()).join(","))).join(","));
    println!("{l}");
}

// cargo test -r -- --exact k9_46 --nocapture --include-ignored
#[test]
#[ignore]
fn k9_46() { 
    let l = InvLink::sinv_knot_from_code(
        [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // good
    assert_eq!(ssi, (0, 2));
}

// cargo test -r -- --exact k15n_103488 --nocapture --include-ignored
#[test]
#[ignore]
fn k15n_103488() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code(
        [[1,11,2,10],[2,20,3,19],[5,17,6,16],[6,25,7,26],[9,22,10,23],[12,30,13,29],[14,8,15,7],[15,27,16,26],[18,4,19,3],[20,11,21,12],[21,1,22,30],[23,4,24,5],[24,18,25,17],[27,8,28,9],[28,14,29,13]],
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // good
    assert_eq!(ssi, (0, 2));
}

// cargo test -r -- --exact knotJ --nocapture --include-ignored
#[test]
#[ignore]
fn knotJ() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code([
        [1,27,2,26],[19,2,20,3],[3,13,4,12],[4,31,5,32],[30,5,31,6],
        [13,7,14,6],[8,27,9,28],[9,1,10,34],[10,18,11,17],[24,11,25,12],
        [14,21,15,22],[28,16,29,15],[33,16,34,17],[18,26,19,25],[20,8,21,7],
        [29,23,30,22],[23,33,24,32]
    ]);

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // good
    assert_eq!(ssi, (0, 2));
}

#[test]
#[ignore]
fn k8_8_and_10_129() { 
    // 8_8a
    let l1 = InvLink::sinv_knot_from_code(
        [[2,8,3,7],[5,15,6,14],[6,4,7,3],[8,11,9,12],[10,1,11,2],[12,18,13,17],[15,5,16,4],[16,14,17,13],[18,9,1,10]],
    );

    // 8_8b
    let l2 = InvLink::sinv_knot_from_code(
        [[1,11,2,10],[2,15,3,16],[4,17,5,18],[5,15,6,14],[7,13,8,12],[9,1,10,18],[11,9,12,8],[13,7,14,6],[16,3,17,4]],
    );

    // 10_129
    let l3 = InvLink::sinv_knot_from_code(
        [[2,15,3,16],[4,17,5,18],[6,19,7,20],[7,13,8,12],[9,15,10,14],[10,2,11,1],[13,9,14,8],[16,5,17,6],[18,3,19,4],[20,12,1,11]],
    );

    let c = P::variable();

    // not good.
    assert_eq!(ssi_invariants(&l1, &c), (0, 0));
    assert_eq!(ssi_invariants(&l2, &c), (0, 0));
    assert_eq!(ssi_invariants(&l3, &c), (0, 0));
}

#[test]
#[ignore]
fn k10_71() { 
    let l = InvLink::sinv_knot_from_code(
        [[1,13,2,12],[3,6,4,7],[5,19,6,18],[7,10,8,11],[8,15,9,16],[11,1,12,22],[13,16,14,17],[14,9,15,10],[17,20,18,21],[19,5,20,4],[21,3,22,2]],
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // not good.
    assert_eq!(ssi, (0, 0));
}

#[test]
#[ignore]
fn k10_104() { 
    let l = InvLink::sinv_knot_from_code(
        [[1,14,2,15],[3,11,4,10],[5,16,6,17],[7,20,8,1],[9,3,10,2],[11,19,12,18],[13,9,14,8],[15,4,16,5],[17,6,18,7],[19,13,20,12]],
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // not good.
    assert_eq!(ssi, (0, 0));

    // khi is not symmetric.
}

#[test]
#[ignore]
fn k13n_1496() { 
    init_logger();
    
    // see [Boyle-Issa,2021]
    let l = InvLink::sinv_knot_from_code(
        [[1,9,2,8],[2,13,3,14],[4,25,5,26],[6,4,7,3],[7,15,8,14],[11,19,12,18],[15,23,16,22],[16,27,17,28],[17,11,18,10],[19,13,20,12],[20,9,21,10],[21,1,22,28],[24,5,25,6],[26,24,27,23]],
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // not good.
    assert_eq!(ssi, (2, 2));
}

#[test]
#[ignore]
fn W4_1() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code(
        [
            [1,9,2,8],[2,13,3,14],[4,25,5,26],[6,4,7,3],[7,15,8,14],
            [11,19,12,18],[15,23,16,22],[16,27,17,28],[17,11,18,10],[19,13,20,12],
            [20,9,21,10],[21,1,22,28],[24,5,25,6],[26,24,27,23]
        ],
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // not good.
    assert_eq!(ssi, (0, 0));
}

#[test]
#[ignore]
// cargo test -r -- --exact W3_1 --nocapture --include-ignored
fn W3_1() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code(
        [
            [3,25,4,24],[4,37,5,38],[7,14,8,15],[9,12,10,13],[11,31,12,30],
            [13,8,14,9],[17,39,18,38],[18,23,19,24],[21,7,22,6],[22,15,23,16],
            [25,3,26,2],[26,19,27,20],[27,34,28,35],[29,32,30,33],[31,11,32,10],
            [33,28,34,29],[35,21,36,20],[36,1,37,2],[39,17,40,16],[40,5,1,6]
        ],
    ).mirror();

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    println!("{ssi:?}");
}

#[test]
#[ignore]
// cargo test -r -- --exact cable4_1 --nocapture --include-ignored
fn cable4_1() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code(
        [
            [1,24,2,25],[2,7,3,8],[5,31,6,30],[6,14,7,13],[9,27,10,26],
            [10,17,11,18],[11,34,12,1],[14,4,15,3],[15,21,16,20],[18,25,19,26],
            [19,8,20,9],[22,30,23,29],[23,13,24,12],[27,16,28,17],[28,33,29,34],
            [31,5,32,4],[32,22,33,21],
        ],
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    println!("{ssi:?}");
}
