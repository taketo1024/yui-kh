use clap::{Parser, Subcommand, ValueEnum};
use log::{info, error};
use simplelog::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    #[arg(long, default_value_t = false)]
    debug: bool
}

#[derive(Subcommand, Debug)]
enum Commands {
    SS {
        name: String,

        #[arg(short, long)]
        link: Option<String>,

        #[arg(short, long)]
        c_value: String,

        #[arg(short = 't', long, default_value_t = ss::CType::i64)]
        c_type: ss::CType,

        #[arg(short, long)]
        output: Option<String>,
    },

    SSBatch {
        targets: String,
        c_value: String,
        
        #[arg(short = 't', long, default_value_t = ss::CType::i64)]
        c_type: ss::CType,

        #[arg(short, long)]
        output: Option<String>,
    },
}

fn main() {
    let args = Cli::parse();

    if args.debug { 
        init_logger();
    }

    let (res, time) = measure(|| 
        guard_panic(||
            match args.command { 
            Commands::SS { name, link, c_value, c_type, output } 
                => ss::run_single(name, link, c_value, c_type, output),
            Commands::SSBatch { targets, c_value, c_type, output }
                => ss::run_batch(targets, c_value, c_type, output)
            }
        )
    );

    if let Ok(res) = res { 
        info!("time: {:?}", time);
        println!("{res}");
        std::process::exit(0)

    } else if let Err(e) = res { 
        error!("{}", e);
        eprintln!("{e}");
        std::process::exit(1)
    }
}

mod ss { 
    use std::str::FromStr;

    use super::*;
    use yui::links::Link;
    use yui::khovanov::invariants::ss::ss_invariant;
    use derive_more::Display;
    use yui::math::traits::{EucRing, EucRingOps};

    #[allow(non_camel_case_types)]
    #[derive(Clone, ValueEnum, Debug, Display)]
    pub enum CType { 
        i64, i128, bigint, 
        gauss, gauss_128, gauss_big,
        eisen, eisen_128, eisen_big,
        none
    }

    pub fn run_single(name: String, link: Option<String>, c_value: String, c_type: CType, output: Option<String>) -> Result<String, Box<dyn std::error::Error>> {
        let l = load_link(&name, &link)?;
        run(&name, &l, &c_value, &c_type, &output)
    }

    pub fn run_batch(targets: String, c_value: String, c_type: CType, output: Option<String>) -> Result<String, Box<dyn std::error::Error>> {
        let data = load_data(&targets)?;
        let mut all_res = String::from("");

        for (name, code) in data { 
            let l = Link::from(&code);
            let res = run(&name, &l, &c_value, &c_type, &output)?;

            if !all_res.is_empty() { 
                all_res += "\n";
            }
            all_res += &format!("{name}: {res}");
        }
        
        Ok(all_res)
    }

    fn run(name: &String, link: &Link, c_value: &String, c_type: &CType, output: &Option<String>) -> Result<String, Box<dyn std::error::Error>> {
        info!("compute ss: {name}, c = {c_value}");

        let res = guard_panic(|| 
            compute_ss_switch(link, c_value, c_type)
        );

        if let Ok(s) = res { 
            info!("{name}: s = {s} (c = {c_value})")
        }

        if let Some(output) = output { 
            write_res(&name, &c_value, &res, &output)?;
        }

        res.map(|s| s.to_string() )
    }

    fn compute_ss_switch(l: &Link, c_value: &String, c_type: &CType) -> Result<i32, Box<dyn std::error::Error>> {
        use num_bigint::BigInt;
        use yui::math::types::quad_int::{GaussInt, EisenInt};

        match c_type { 
            CType::i64        => compute_ss::<i64>(l, c_value),
            CType::i128       => compute_ss::<i128>(l, c_value),
            CType::bigint     => compute_ss::<BigInt>(l, c_value),
            CType::gauss      => compute_ss::<GaussInt<i64>>(l, c_value), 
            CType::gauss_128  => compute_ss::<GaussInt<i128>>(l, c_value), 
            CType::gauss_big  => compute_ss::<GaussInt<BigInt>>(l, c_value), 
            CType::eisen      => compute_ss::<EisenInt<i64>>(l, c_value), 
            CType::eisen_128  => compute_ss::<EisenInt<i128>>(l, c_value), 
            CType::eisen_big  => compute_ss::<EisenInt<BigInt>>(l, c_value), 
            _ => err!("")
        }
    }

    fn compute_ss<R>(l: &Link, c_value: &String) -> Result<i32, Box<dyn std::error::Error>>
    where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
        if let Ok(c) = c_value.parse::<R>() { 
            Ok( ss_invariant(l, &c, true) )
        } else { 
            err!("cannot parse c: '{}' as type: {}.", c_value, std::any::type_name::<R>())
        }
    }

    fn write_res(name: &String, c_value: &String, res: &Result<i32, Box<dyn std::error::Error>>, output: &String) -> Result<(), Box<dyn std::error::Error>> { 
        use std::fs::OpenOptions;
        use std::path::Path;

        let file_exists = Path::new(output).exists();
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .append(true)
            .open(output)?;
        let mut wtr = csv::Writer::from_writer(file);

        if !file_exists { 
            wtr.write_record(vec!["name", "c", "ss"])?;
        }

        let res = if let Ok(s) = res { 
            s.to_string()
        } else { 
            "!".to_string()
        };

        wtr.write_record(vec![name, c_value, &res])?;
        wtr.flush()?;

        info!("write: {output}");

        Ok(())
    }
}

use indexmap::IndexMap;
use yui::links::{Link, links::Edge};

type PDCode = Vec<[Edge; 4]>;
type Data = IndexMap<String, PDCode>;

fn load_data(path: &String) -> Result<Data, Box<dyn std::error::Error>> { 
    let json = std::fs::read_to_string(path)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}

fn load_link(name: &String, pd_code: &Option<String>) -> Result<Link, Box<dyn std::error::Error>> { 
    if let Some(pd_code) = pd_code { 
        let pd_code: PDCode = serde_json::from_str(&pd_code)?;
        let l = Link::from(&pd_code);
        Ok(l)
    } else if let Ok(l) = Link::load(&name) { 
        Ok(l)
    } else {
        err!("cannot load {name}")
    }
}

fn init_logger() {
    TermLogger::init(
        LevelFilter::Info,
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto
    ).unwrap()
}

fn measure<F, Res>(proc: F) -> (Res, std::time::Duration) 
where F: FnOnce() -> Res { 
    let start = std::time::Instant::now();
    let res = proc();
    let time = start.elapsed();
    (res, time)
}

fn guard_panic<F, R>(f: F) -> Result<R, Box<dyn std::error::Error>>
where F: FnOnce() -> Result<R, Box<dyn std::error::Error>> + std::panic::UnwindSafe {
    std::panic::catch_unwind(|| {
        f()
    }).unwrap_or_else(|e| {
        let info = match e.downcast::<String>() {
            Ok(v) => *v,
            Err(e) => match e.downcast::<&str>() {
                Ok(v) => v.to_string(),
                _ => "Unknown Source of Error".to_owned()
            }
        };
        err!("panic: {info}")
    })
}

#[derive(Debug, derive_more::Display)]
struct Error(String);
impl std::error::Error for Error {}

macro_rules! err {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        Err( Error(msg).into() )
    }}
}
pub(crate) use err;