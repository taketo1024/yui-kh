use log::{info, error};
use clap::{Parser, Subcommand};
use super::cmd::{kh, ckh, ss};
use super::utils::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct CliArgs {
    #[command(subcommand)]
    pub command: Cmd
}

#[derive(Subcommand, Debug)]
pub enum Cmd {
    Kh(kh::Args),
    Ckh(ckh::Args),
    SS(ss::Args),
}

impl CliArgs { 
    fn log_level(&self) -> log::LevelFilter { 
        use log::LevelFilter::*;
        let level = match &self.command { 
            Cmd::Kh(args)  => args.log,
            Cmd::Ckh(args) => args.log,
            Cmd::SS(args)  => args.log,
        };
        match level {
            1 => Info,
            2 => Debug,
            3 => Trace,
            _ => Off,
        }
    }
}

pub struct App {
    pub args: CliArgs
}

impl App { 
    pub fn new() -> Self { 
        let args = CliArgs::parse();
        App { args }
    }

    pub fn run(&self) -> Result<String, i32> { 
        self.init_logger();

        info!("args: {:?}", self.args);
        info!("int-type: {}", std::any::type_name::<super::utils::dispatch::Int>());

        let (res, time) = measure(||
            self.dispatch()
        );

        let res = res.map_err(|e| { 
            error!("{}", e);
            eprintln!("\x1b[0;31merror\x1b[0m: {e}");
            1 // error code
        });

        info!("time: {:?}", time);

        res
    }

    fn init_logger(&self) {
        use simplelog::*;
        TermLogger::init(
            self.args.log_level(),
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto
        ).unwrap()
    }

    fn dispatch(&self) -> Result<String, Box<dyn std::error::Error>> { 
        guard_panic(||
            match &self.args.command { 
                Cmd::Kh(args)  => kh::run(args),
                Cmd::Ckh(args) => ckh::run(args),
                Cmd::SS(args)  => ss::run(args)
            }
        )
    }
}