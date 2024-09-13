use std::marker::PhantomData;
use std::str::FromStr;
use yui::{EucRing, EucRingOps};
use yui_homology::{DisplayForGrid, DisplaySeq, DisplayTable, GridDeg, GridTrait, RModStr, XHomologyBase};
use yui_kh::{KhComplex, KhGen};
use crate::app::utils::*;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(App, args)
}

#[derive(Clone, Default, Debug, clap::Args)]
pub struct Args { 
    pub link: String,

    #[arg(short = 't', long, default_value = "Z")]
    pub c_type: CType,

    #[arg(short, long, default_value = "0")]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub reduced: bool,

    #[arg(short = 'g', long)]
    pub show_generators: bool,

    #[arg(long)]
    pub no_simplify: bool,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

pub struct App<R>
where
    R: EucRing + FromStr,
    for<'x> &'x R: EucRingOps<R>,
{
    args: Args,
    _ring: PhantomData<R>
}

impl<R> App<R>
where
    R: EucRing + FromStr,
    for<'x> &'x R: EucRingOps<R>,
{
    pub fn new(args: Args) -> Self { 
        App { args, _ring: PhantomData }
    }

    pub fn run(&self) -> Result<String, Box<dyn std::error::Error>> { 
        let args = &self.args;
        let (h, t) = parse_pair::<R>(&args.c_value)?;
    
        let bigraded = h.is_zero() && t.is_zero();
        let poly = ["H", "0,T"].contains(&args.c_value.as_str());
        let show_gens = args.show_generators;
    
        if args.reduced && !t.is_zero() { 
            return err!("{t} != 0 is not allowed for reduced.");
        }
    
        let l = load_link(&args.link, args.mirror)?;
        let ckh = if args.no_simplify { 
            KhComplex::new_v1(&l, &h, &t, args.reduced)
        } else {
            KhComplex::new(&l, &h, &t, args.reduced)
        };
    
        let mut b = string_builder::Builder::new(1024);
    
        if bigraded { 
            let ckh = ckh.into_bigraded();
            let kh = ckh.homology(show_gens);
    
            b.append(kh.display_table("i", "j"));
    
            if show_gens { 
                b.append(display_gens(kh.inner()));
            }
        } else if poly { 
            let kh = ckh.homology(true);
            let gens = kh.gen_table();
    
            b.append(gens.display_table("i", "j"));
    
            if show_gens { 
                b.append(display_gens(kh.inner()));
            }
        } else { 
            let kh = ckh.homology(show_gens);
    
            b.append(kh.display_seq("i"));
    
            if show_gens { 
                b.append(display_gens(kh.inner()));
            }
        };
    
        let res = b.string()?;
        Ok(res)
    }
}

fn display_gens<I, R>(kh: &XHomologyBase<I, KhGen, R>) -> String
where I: GridDeg, R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let mut res = String::new();

    for i in kh.support() {
        let h = &kh[i];
        if h.is_zero() { continue }

        res += "\n";
        res += &format!("{i}: {}\n", h.display_for_grid());

        let r = h.rank() + h.tors().len();
        for i in 0..r { 
            let z = h.gen_chain(i);
            res += &format!("  {i}: {z}\n");
        }
    }

    res
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
        let res = dispatch(&args);
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
        let res = dispatch(&args);
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
            let res = dispatch(&args);
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
            let res = dispatch(&args);
            assert!(res.is_ok());
        }
    }
}