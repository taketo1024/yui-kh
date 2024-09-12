use std::str::FromStr;
use yui::{EucRing, EucRingOps};
use yui_homology::{DisplayTable, DisplaySeq};
use yui_kh::KhComplex;
use crate::app::utils::*;

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

    #[arg(long)]
    no_simplify: bool,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

pub fn run(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(&args.c_value, &args.c_type, compute_homology, args)
}

fn compute_homology<R>(args: &Args) -> Result<String, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let (h, t) = parse_pair::<R>(&args.c_value)?;
    let bigraded = h.is_zero() && t.is_zero();

    if args.reduced && !t.is_zero() { 
        return err!("{t} != 0 is not allowed for reduced.");
    }

    let l = load_link(&args.link, args.mirror)?;
    let ckh = if args.no_simplify { 
        KhComplex::new_v1(&l, &h, &t, args.reduced)
    } else {
        KhComplex::new(&l, &h, &t, args.reduced)
    };

    let res = if bigraded { 
        let ckh = ckh.into_bigraded();
        let kh = ckh.homology(false);
        kh.display_table("i", "j")
    } else { 
        let kh = ckh.homology(false);
        kh.display_seq("i")
    };

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
            c_value: "0".to_string(),
            c_type: CType::Z,
            mirror: true,
            reduced: true,
            ..Default::default()
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[cfg(feature = "poly")]
    mod poly_tests { 
        use super::*;
        
        #[test]
        fn test_qpoly_h() { 
            let args = Args {
                link: "3_1".to_string(),
                c_value: "H".to_string(),
                c_type: CType::Q,
                ..Default::default()
            };
            let res = run(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_qpoly_t() { 
            let args = Args {
                link: "3_1".to_string(),
                c_value: "0,T".to_string(),
                c_type: CType::Q,
                ..Default::default()
            };
            let res = run(&args);
            assert!(res.is_ok());
        }
    }
}