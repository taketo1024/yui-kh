#![allow(non_snake_case)]
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

// cargo test -r -- --exact k9_46 --nocapture --include-ignored
#[test]
#[ignore]
fn k9_46() { 
    let l = InvLink::from_code(
        [[17,7,0,6],[12,5,13,6],[11,1,12,0],[7,17,8,16],[4,13,5,14],[1,11,2,10],[15,9,16,8],[14,3,15,4],[9,3,10,2]], 
        [(6,12),(7,11),(1,17),(5,13),(2,16),(8,10),(4,14),(3,15)],
        Some(0)
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
    
    let l = InvLink::from_code(
        [[1,11,2,10],[2,20,3,19],[5,17,6,16],[6,25,7,26],[9,22,10,23],[12,30,13,29],[14,8,15,7],[15,27,16,26],[18,4,19,3],[20,11,21,12],[21,1,22,30],[23,4,24,5],[24,18,25,17],[27,8,28,9],[28,14,29,13]],
        [(2,30),(3,29),(4,28),(5,27),(6,26),(7,25),(8,24),(9,23),(10,22),(11,21),(12,20),(13,19),(14,18),(15,17)],
        Some(1)
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
    
    let l = InvLink::from_code(
        [
            [0,26,1,25],[18,1,19,2],[2,12,3,11],[3,30,4,31],
            [29,4,30,5],[12,6,13,5],[7,26,8,27],[8,0,9,33],
            [9,17,10,16],[23,10,24,11],[13,20,14,21],[27,15,28,14],
            [32,15,33,16],[17,25,18,24],[19,7,20,6],[28,22,29,21],
            [22,32,23,31]
        ], 
        [(1,33),(2,32),(3,31),(4,30),(5,29),(6,28),(7,27),(8,26),(9,25),(10,24),(11,23),(12,22),(13,21),(14,20),(15,19),(16,18)],
        Some(0)
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // good
    assert_eq!(ssi, (0, 2));
}

#[test]
#[ignore]
fn k8_8_and_10_129() { 
    // 8_8a
    let l1 = InvLink::from_code(
        [[2,8,3,7],[5,15,6,14],[6,4,7,3],[8,11,9,12],[10,1,11,2],[12,18,13,17],[15,5,16,4],[16,14,17,13],[18,9,1,10]],
        [(2,18),(3,17),(4,16),(5,15),(6,14),(7,13),(8,12),(9,11)],
        Some(1)
    );

    // 8_8b
    let l2 = InvLink::from_code(
        [[1,11,2,10],[2,15,3,16],[4,17,5,18],[5,15,6,14],[7,13,8,12],[9,1,10,18],[11,9,12,8],[13,7,14,6],[16,3,17,4]],
        [(2,18),(3,17),(4,16),(5,15),(6,14),(7,13),(8,12),(9,11)],
        Some(1)
    );

    // 10_129
    let l3 = InvLink::from_code(
        [[2,15,3,16],[4,17,5,18],[6,19,7,20],[7,13,8,12],[9,15,10,14],[10,2,11,1],[13,9,14,8],[16,5,17,6],[18,3,19,4],[20,12,1,11]],
        [(2,20),(3,19),(4,18),(5,17),(6,16),(7,15),(8,14),(9,13),(10,12)],
        Some(1)
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
    let l = InvLink::from_code(
        [[1,13,2,12],[3,6,4,7],[5,19,6,18],[7,10,8,11],[8,15,9,16],[11,1,12,22],[13,16,14,17],[14,9,15,10],[17,20,18,21],[19,5,20,4],[21,3,22,2]],
        [(2,22),(3,21),(4,20),(5,19),(6,18),(7,17),(8,16),(9,15),(10,14),(11,13)],
        Some(1)
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // not good.
    assert_eq!(ssi, (0, 0));
}

#[test]
#[ignore]
fn k10_104() { 
    let l = InvLink::from_code(
        [[1,14,2,15],[3,11,4,10],[5,16,6,17],[7,20,8,1],[9,3,10,2],[11,19,12,18],[13,9,14,8],[15,4,16,5],[17,6,18,7],[19,13,20,12]],
        [(2,20),(3,19),(4,18),(5,17),(6,16),(7,15),(8,14),(9,13),(10,12)],
        Some(1)
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
    let l = InvLink::from_code(
        [[1,9,2,8],[2,13,3,14],[4,25,5,26],[6,4,7,3],[7,15,8,14],[11,19,12,18],[15,23,16,22],[16,27,17,28],[17,11,18,10],[19,13,20,12],[20,9,21,10],[21,1,22,28],[24,5,25,6],[26,24,27,23]],
        [(2,28),(3,27),(4,26),(5,25),(6,24),(7,23),(8,22),(9,21),(10,20),(11,19),(12,18),(13,17),(14,16)],
        Some(1)
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
    
    let l = InvLink::from_code(
        [[1,9,2,8],[2,13,3,14],[4,25,5,26],[6,4,7,3],[7,15,8,14],[11,19,12,18],[15,23,16,22],[16,27,17,28],[17,11,18,10],[19,13,20,12],[20,9,21,10],[21,1,22,28],[24,5,25,6],[26,24,27,23]],
        [(2,28),(3,27),(4,26),(5,25),(6,24),(7,23),(8,22),(9,21),(10,20),(11,19),(12,18),(13,17),(14,16)],
        Some(1)
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    // not good.
    assert_eq!(ssi, (0, 0));
}

#[test]
#[ignore]
// cargo test -r -- --exact knotQ --nocapture --include-ignored
fn knotQ() { 
    init_logger();
    
    let l = InvLink::from_code(
        [
            [1,11,2,10],[5,27,6,26],[6,35,7,36],[7,14,8,15],[11,3,12,2],
            [12,19,13,20],[13,8,14,9],[15,37,16,36],[16,25,17,26],[20,9,21,10],
            [21,33,22,32],[23,5,24,4],[24,17,25,18],[27,34,28,35],[29,23,30,22],
            [30,39,31,40],[33,28,34,29],[37,19,38,18],[38,3,39,4],[40,31,1,32]
        ],
        [
            (2,40),(3,39),(4,38),(5,37),(6,36),
            (7,35),(8,34),(9,33),(10,32),(11,31),
            (12,30),(13,29),(14,28),(15,27),(16,26),
            (17,25),(18,24),(19,23),(18,22),(19,21)
        ],
        Some(1)
    );

    let c = P::variable();
    let ssi = ssi_invariants(&l, &c);

    println!("{ssi:?}");
}
