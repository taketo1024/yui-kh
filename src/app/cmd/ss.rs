use std::str::FromStr;
use indexmap::IndexMap;
use log::{info, error};
use yui_link::{Link, Edge};
use yui_core::{EucRing, EucRingOps};
use yui_khovanov::ss_invariant;
use crate::utils::*;

#[derive(Debug, clap::Args)]
pub struct Args { 
    name: String,

    #[arg(short, long)]
    link: Option<String>,

    #[arg(short, long)]
    c_value: String,

    #[arg(short = 't', long, default_value = "z")]
    c_type: CType,
}

#[derive(Debug, clap::Args)]
pub struct BatchArgs { 
    #[arg(short, long)]
    c_value: String,
    
    #[arg(short = 't', long, default_value = "z")]
    c_type: CType,

    #[arg(short, long)]
    data: String,

    #[arg(short, long)]
    output: Option<String>
}

pub fn run(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    info!("compute ss: {}, c = {}", args.name, args.c_value);

    let l = load_link(&args.name, &args.link, false)?;
    let s = guard_panic(|| 
        dispatch_eucring!(&args.c_type, compute_ss, &l, &args.c_value)
    )?;

    info!("{}: ss = {} (c = {})", args.name, s, args.c_value);

    Ok(s.to_string())
}

pub fn run_batch(args: &BatchArgs) -> Result<String, Box<dyn std::error::Error>> {
    let data: IndexMap<String, Vec<[Edge; 4]>> = load_json(&args.data)?;
    let mut all_res = String::from("");

    for (name, code) in data { 
        let l = Link::from(&code);
        let res = guard_panic(|| 
            dispatch_eucring!(&args.c_type, compute_ss, &l, &args.c_value)
        );

        let s = if let Ok(s) = res { 
            s.to_string()
        } else { 
            if let Err(e) = res {
                error!("{}", e);
                eprintln!("{e}");
            }
            "!".to_string()            
        };
        
        if let Some(output) = &args.output { 
            write_csv(output, vec![&name, &s])?;
        }

        if !all_res.is_empty() { 
            all_res += "\n";
        }
        
        all_res += &format!("{name}: {s}");
    }
    
    Ok(all_res)
}

fn compute_ss<R>(l: &Link, c_value: &String) -> Result<i32, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let Ok(c) = R::from_str(c_value) else { 
        return err!("cannot parse c: '{}' as type: {}.", c_value, std::any::type_name::<R>())
    };

    let ss = ss_invariant(l, &c, true);
    Ok(ss)
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let args = Args {
        	name: "3_1".to_string(),
            link: None,
        	c_value: "2".to_string(),
        	c_type: CType::Z,
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let args = Args {
        	name: "4_1".to_string(),
            link: Some("[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string()),
        	c_value: "3".to_string(),
        	c_type: CType::Z,
        };
        let res = run(&args);
        assert!(res.is_ok());
    }
}