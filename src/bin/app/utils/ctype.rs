use clap::ValueEnum;
use derive_more::Display;

#[derive(Clone, ValueEnum, Display, Debug, Default)]
#[clap(rename_all="verbatim")]
pub enum CType { 
    #[default] Z, 
    Q, F2, F3, Gauss, Eisen
}