use std::marker::PhantomData;
use std::str::FromStr;
use log::info;
use yui::{EucRing, EucRingOps};
use yui_kh::{ss_invariant, ss_invariant_v1};
use crate::app::utils::*;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(App, args)
}

#[derive(Clone, Default, Debug, clap::Args)]
pub struct Args { 
    pub link: String,

    #[arg(short = 't', long, default_value = "Z")]
    pub c_type: CType,

    #[arg(short, long)]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub unreduced: bool,

    #[arg(long)]
    pub old: bool,

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

        info!("compute ss: c = {}", args.c_value);

        let l = load_link(&args.link, args.mirror)?;
        let Ok(c) = R::from_str(&args.c_value) else { 
            return err!("cannot parse c: '{}' as type: {}.", &args.c_value, std::any::type_name::<R>())
        };
    
        let ss = if args.old { 
            ss_invariant_v1(&l, &c, !args.unreduced)
        } else { 
            ss_invariant(&l, &c, !args.unreduced)
        };
    
        Ok(ss.to_string())
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let args = Args {
            link: "3_1".to_string(),
        	c_value: "2".to_string(),
            ..Default::default()
        };
        let res = dispatch(&args);

        assert!(res.is_ok());
        assert_eq!(res.unwrap(), "-2".to_string());
    }

    #[test]
    fn test2() { 
        let args = Args {
        	link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
        	c_value: "3".to_string(),
            ..Default::default()
        };
        let res = dispatch(&args);

        assert!(res.is_ok());
        assert_eq!(res.unwrap(), "-2".to_string());
    }

    #[cfg(feature = "poly")]
    mod poly_tests { 
        use super::*;
        #[test]
        fn test1() { 
            let args = Args {
                link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
                c_value: "H".to_string(),
                c_type: CType::Q,
                ..Default::default()
            };
            let res = dispatch(&args);
            
            assert!(res.is_ok());
            assert_eq!(res.unwrap(), "-2".to_string());
        }    
    }
}