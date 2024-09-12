use crate::app::utils::*;
use itertools::Itertools;
use std::str::FromStr;
use yui::{Ring, RingOps};
use yui_homology::{ChainComplexCommon, DisplaySeq, DisplayTable};
use yui_kh::KhComplex;

#[derive(Debug, clap::Args, Default)]
pub struct Args {
    link: String,

    #[arg(short, long, default_value = "0")]
    c_value: String,

    #[arg(short = 't', long, default_value = "Z")]
    c_type: CType,

    #[arg(short, long)]
    mirror: bool,

    #[arg(short, long)]
    reduced: bool,

    #[arg(short = 'a', long)]
    with_alpha: bool,

    #[arg(long)]
    no_simplify: bool,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

pub fn run(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_ring!(&args.c_value, &args.c_type, describe_ckh, args)
}

fn describe_ckh<R>(args: &Args) -> Result<String, Box<dyn std::error::Error>>
where
    R: Ring + FromStr,
    for<'x> &'x R: RingOps<R>,
{
    let (h, t) = parse_pair::<R>(&args.c_value)?;
    let poly = ["H", "0,T", "H,T"].contains(&args.c_value.as_str());
    let bigraded = h.is_zero() && t.is_zero();

    if args.reduced && !t.is_zero() {
        return err!("{t} != 0 is not allowed for reduced.");
    }

    let l = load_link(&args.link, args.mirror)?;
    let c = if args.no_simplify {
        KhComplex::new_v1(&l, &h, &t, args.reduced)
    } else {
        KhComplex::new(&l, &h, &t, args.reduced)
    };

    let vs = if args.with_alpha {
        c.canon_cycles()
            .iter()
            .map(|z| (0, c[0].vectorize(z)))
            .collect_vec()
    } else {
        vec![]
    };

    let mut b = string_builder::Builder::new(1024);

    if bigraded {
        let c = c.into_bigraded();
        b.append(c.display_table("i", "j") + "\n");
        b.append(c.display_d() + "\n");
    } else if poly {
        b.append(c.gen_table().display_table("i", "j") + "\n");
        b.append(c.display_d() + "\n");
    } else {
        b.append(c.display_seq("i") + "\n");
        b.append(c.display_d() + "\n");
    }

    if args.with_alpha {
        b.append("\n");
        for (i, (_, v)) in vs.iter().enumerate() {
            b.append(format!(
                "a[{i}] = [{}]\n",
                v.to_dense()
                    .iter()
                    .map(|r| r.to_string())
                    .collect_vec()
                    .join(", ")
            ));
        }
    }

    let res = b.string()?;
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let args = Args {
            link: "3_1".to_string(),
            c_value: "0".to_string(),
            ..Default::default()
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() {
        let args = Args {
            link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
            c_value: "2".to_string(),
            mirror: true,
            reduced: true,
            with_alpha: true,
            ..Default::default()
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[cfg(feature = "poly")]
    mod poly_tests {
        use super::*;

        #[test]
        fn test_zpoly_h() {
            let args = Args {
                link: "3_1".to_string(),
                c_value: "H".to_string(),
                c_type: CType::Z,
                ..Default::default()
            };
            let res = run(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_zpoly_t() {
            let args = Args {
                link: "3_1".to_string(),
                c_value: "T".to_string(),
                c_type: CType::Z,
                ..Default::default()
            };
            let res = run(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_zpoly_ht() {
            let args = Args {
                link: "3_1".to_string(),
                c_value: "H,T".to_string(),
                c_type: CType::Z,
                ..Default::default()
            };
            let res = run(&args);
            assert!(res.is_ok());
        }
    }
}
