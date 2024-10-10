#![allow(non_snake_case, unused)]
use itertools::Itertools;
use num_traits::Zero;

use log::*;
use yui::poly::HPoly;
use yui::FF2;
use yui_homology::DisplayTable;
use yui_link::InvLink;
use yui_kh::khi::{ssi_invariants, KhIComplex};

type R = FF2;
type P = HPoly<'H', R>;

fn init_logger(l: LevelFilter) { 
    use simplelog::*;

    let mut cb = simplelog::ConfigBuilder::new();
    cb.set_location_level(LevelFilter::Off);
    cb.set_target_level(LevelFilter::Off);
    cb.set_thread_level(LevelFilter::Off);
    cb.set_level_color(Level::Trace, Some(Color::Green));
    let config = cb.build();

    TermLogger::init(
        l,
        config,
        TerminalMode::Mixed,
        ColorChoice::Always
    ).unwrap();
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
    init_logger(LevelFilter::Info);
    
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
    init_logger(LevelFilter::Info);
    
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
    init_logger(LevelFilter::Info);
    
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
    init_logger(LevelFilter::Info);
    
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
    init_logger(LevelFilter::Info);
    
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
    init_logger(LevelFilter::Info);
    
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

#[test]
#[ignore]
// cargo test -r -- --exact k9_46_interlock --nocapture --include-ignored
fn k9_46_interlock() { 
    type R = HPoly<'H', FF2>;

    init_logger(log::LevelFilter::Info);

    let l = InvLink::sinv_knot_from_code([
        [2,50,3,49],[4,48,5,47],[7,24,8,25],[9,52,10,53],[12,60,13,59],
        [14,58,15,57],[16,36,17,35],[18,34,19,33],[20,32,21,31],[21,40,22,41],
        [23,38,24,39],[25,6,26,7],[26,46,27,45],[28,44,29,43],[30,42,31,41],
        [32,20,33,19],[34,18,35,17],[37,54,38,55],[39,22,40,23],[42,30,43,29],
        [44,28,45,27],[46,6,47,5],[48,4,49,3],[50,2,51,1],[51,10,52,11],
        [53,8,54,9],[55,36,56,37],[56,16,57,15],[58,14,59,13],[60,12,1,11]
    ]);
    let (h, _t) = (R::variable(), R::zero());
    let c = KhIComplex::new_v2(&l, &h, false);
    let h = c.homology();

    h.gen_table().print_table("i", "j");
}

#[test]
#[ignore]
// cargo test -r -- --exact k9_46_interlock3 --nocapture --include-ignored
fn k9_46_interlock3() { 
    type R = HPoly<'H', FF2>;

    init_logger(log::LevelFilter::Info);

    let l = InvLink::sinv_knot_from_code([
        [2,84,3,83],[4,82,5,81],[7,26,8,27],[9,86,10,87],[12,94,13,93],
        [14,92,15,91],[16,68,17,67],[18,66,19,65],[20,64,21,63],[21,42,22,43],
        [23,72,24,73],[25,70,26,71],[27,6,28,7],[28,80,29,79],[30,78,31,77],
        [32,76,33,75],[34,52,35,51],[36,50,37,49],[39,56,40,57],[41,22,42,23],
        [44,62,45,61],[46,60,47,59],[48,38,49,37],[50,36,51,35],[52,34,53,33],
        [53,74,54,75],[55,40,56,41],[57,38,58,39],[58,48,59,47],[60,46,61,45],
        [62,44,63,43],[64,20,65,19],[66,18,67,17],[69,88,70,89],[71,24,72,25],
        [73,54,74,55],[76,32,77,31],[78,30,79,29],[80,6,81,5],[82,4,83,3],
        [84,2,85,1],[85,10,86,11],[87,8,88,9],[89,68,90,69],[90,16,91,15],
        [92,14,93,13],[94,12,1,11]
    ]);
    let (h, _t) = (R::variable(), R::zero());
    let c = KhIComplex::new_v2(&l, &h, false);
    let h = c.homology();

    h.gen_table().print_table("i", "j");
}

#[test]
#[ignore]
// cargo test -r -- --exact knotJ_interlock --nocapture --include-ignored
fn knotJ_interlock() { 
    type R = HPoly<'H', FF2>;

    init_logger(log::LevelFilter::Info);

    let l = InvLink::sinv_knot_from_code([
        [1,17,2,16],[6,23,7,24],[7,115,8,114],[8,103,9,104],[9,3,10,2],
        [10,96,11,95],[14,47,15,48],[15,91,16,90],[18,113,19,114],[19,25,20,24],
        [20,5,21,6],[21,101,22,100],[26,112,27,111],[28,33,29,34],[29,73,30,72],
        [31,107,32,106],[34,72,35,71],[36,63,37,64],[37,43,38,42],[38,55,39,56],
        [39,83,40,82],[44,77,45,78],[45,61,46,60],[48,13,49,14],[49,93,50,92],
        [50,88,51,87],[56,41,57,42],[57,65,58,64],[58,85,59,86],[59,53,60,52],
        [61,77,62,76],[66,83,67,84],[67,55,68,54],[68,43,69,44],[69,63,70,62],
        [70,36,71,35],[74,107,75,108],[75,31,76,30],[78,53,79,54],[79,85,80,84],
        [80,65,81,66],[81,41,82,40],[86,52,87,51],[88,93,89,94],[89,13,90,12],
        [91,47,92,46],[94,12,95,11],[96,3,97,4],[97,103,98,102],[98,115,99,116],
        [99,23,100,22],[104,17,105,18],[105,1,106,120],[108,73,109,74],[109,33,110,32],
        [110,28,111,27],[116,101,117,102],[117,5,118,4],[118,25,119,26],[119,113,120,112]
    ]);
    let (h, _t) = (R::variable(), R::zero());
    let c = KhIComplex::new_v2(&l, &h, false);
    let h = c.homology();

    h.gen_table().print_table("i", "j");
}
