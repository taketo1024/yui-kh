use clap::ValueEnum;
use derive_more::Display;
use serde::Deserialize;

#[derive(Clone, ValueEnum, Display, Debug, Deserialize, Default)]
#[clap(rename_all="verbatim")]
pub enum CType { 
    #[default] Z, 
    Q, F2, F3, Gauss, Eisen
}