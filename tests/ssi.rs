use yui::poly::HPoly;
use yui::FF2;
use yui_kh::khi::ssi_invariants;
use yui_link::InvLink;

macro_rules! test {
    ($(#[$m:meta])* $test:ident, $name:literal, $expected:expr) => {
        $(#[$m])* 
        #[test]
        fn $test() -> Result<(), Box<dyn std::error::Error>> { 
            type R = FF2;
            type P = HPoly<'H', R>;
            let c = P::variable();

            let l = InvLink::load($name)?;
            let ssi = ssi_invariants(&l, &c);
            assert_eq!(ssi, $expected);

            Ok(())
        }
    }
}

test!(k3_1, "3_1", (2, 2));
test!(k4_1, "4_1", (0, 0));
test!(k5_1, "5_1", (4, 4));
test!(k5_2a, "5_2a", (2, 2));
test!(k5_2b, "5_2b", (2, 2));
test!(k6_1a, "6_1a", (0, 0));
test!(k6_1b, "6_1b", (0, 0));
test!(k6_2a, "6_2a", (2, 2));
test!(k6_2b, "6_2b", (2, 2));
test!(k6_3, "6_3", (0, 0));
test!(k7_1, "7_1", (6, 6));
test!(k7_2a, "7_2a", (2, 2));
test!(k7_2b, "7_2b", (2, 2));
test!(k7_3a, "7_3a", (4, 4));
test!(k7_3b, "7_3b", (4, 4));
test!(k7_4a, "7_4a", (2, 2));
test!(k7_4b, "7_4b", (2, 2));
test!(k7_5a, "7_5a", (4, 4));
test!(k7_5b, "7_5b", (4, 4));
test!(k7_6a, "7_6a", (-2, -2));
test!(k7_6b, "7_6b", (-2, -2));
test!(k7_7a, "7_7a", (0, 0));
test!(k7_7b, "7_7b", (0, 0));