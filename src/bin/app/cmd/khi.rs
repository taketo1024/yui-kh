use std::marker::PhantomData;
use std::str::FromStr;
use yui::{EucRing, EucRingOps};
use yui_homology::{DisplayForGrid, DisplaySeq, DisplayTable, GridTrait, RModStr};
use yui_kh::khi::{KhIChain, KhIChainExt, KhIComplex, KhIHomology};
use yui_link::InvLink;
use crate::app::utils::*;
use crate::app::err::*;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(App, args)
}

#[derive(Clone, Default, Debug, clap::Args)]
pub struct Args { 
    pub link: String,

    #[arg(short = 't', long, default_value = "F2")]
    pub c_type: CType,

    #[arg(short, long, default_value = "0")]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub reduced: bool,

    #[arg(short = 'g', long)]
    pub show_gens: bool,

    #[arg(short = 'a', long)]
    pub show_alpha: bool,

    #[arg(short = 's', long)]
    pub show_ssi: bool,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

pub struct App<R>
where
    R: EucRing + FromStr,
    for<'x> &'x R: EucRingOps<R>,
{
    args: Args,
    buff: String,
    _ring: PhantomData<R>
}

impl<R> App<R>
where
    R: EucRing + FromStr,
    for<'x> &'x R: EucRingOps<R>,
{
    pub fn new(args: Args) -> Self { 
        let buff = String::with_capacity(1024);
        App { args, buff, _ring: PhantomData }
    }

    pub fn run(&mut self) -> Result<String, Box<dyn std::error::Error>> { 
        let (h, t) = parse_pair::<R>(&self.args.c_value)?;

        ensure!(self.args.c_type == CType::F2, "Only `-t F2` is supported.");

        if self.args.reduced { 
            ensure!(t.is_zero(), "`t` must be zero for reduced.");
        }
        if self.args.show_alpha { 
            ensure!(t.is_zero(), "`t` must be zero to have alpha.");
        }
        if self.args.show_ssi { 
            ensure!(!h.is_zero() && !h.is_unit(), "`h` must be non-zero, non-invertible to compute ssi.");
            ensure!(t.is_zero(), "`t` must be zero to compute ss.");
        }
    
        let l = load_sinv_knot(&self.args.link, self.args.mirror)?;
        let ckhi = KhIComplex::<R>::new(&l, &h, &t, self.args.reduced);
        let khi = ckhi.homology();

        let bigraded = h.is_zero() && t.is_zero() || 
            ["H", "0,T"].contains(&self.args.c_value.as_str());

        if bigraded { 
            let khi = khi.clone().into_bigraded();
            self.out(&khi.display_table("i", "j"));
        } else { 
            self.out(&khi.display_seq("i"));
        }

        if self.args.show_gens { 
            self.show_gens(&khi);
        }

        if self.args.show_alpha { 
            let zs = ckhi.canon_cycles();
            self.show_alpha(&khi, zs);
        }

        if self.args.show_ssi { 
            let zs = ckhi.canon_cycles();
            self.show_ssi(&l, &h, &khi, zs);
        }

        Ok(self.flush())
    }

    fn show_gens(&mut self, khi: &KhIHomology<R>) { 
        for i in khi.support() {
            let h = &khi[i];
            if h.is_zero() { continue }

            self.out(&format!("KhI[{i}]: {}", h.display_for_grid()));

            let r = h.rank() + h.tors().len();
            for i in 0..r { 
                let z = h.gen_chain(i);
                self.out(&format!("  {i}: {z}"));
            }
            self.out("");
        }
    }

    fn show_alpha(&mut self, khi: &KhIHomology<R>, zs: &[KhIChain<R>]) { 
        for (i, z) in zs.iter().enumerate() { 
            let v = khi[z.h_deg()].vectorize_euc(z);
            self.out(&format!("a[{i}] in KhI[{}]: {}", z.h_deg(), vec2str(&v)));
        }
        self.out("");
    }

    fn show_ssi(&mut self, l: &InvLink, c: &R, khi: &KhIHomology<R>, zs: &[KhIChain<R>]) { 
        assert!(!c.is_unit() && !c.is_unit());

        use yui_kh::misc::div_vec;

        let l = l.link();
        let w = l.writhe();
        let r = l.seifert_circles().len() as i32;

        for (i, z) in zs.iter().enumerate() { 
            let h = &khi[z.h_deg()];
            let v = h.vectorize(z).subvec(0..h.rank());
            let d = div_vec(&v, c).unwrap();
            let s = 2 * d + w - r + 1;

            self.out(&format!("ss[{i}] = {s} (d = {d}, w = {w}, r = {r})"));
        }
    }

    fn out(&mut self, str: &str) { 
        self.buff.push_str(str);
        self.buff.push('\n');
    }

    fn flush(&mut self) -> String { 
        let res = std::mem::take(&mut self.buff);
        res.trim().to_string()
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let args = Args { 
            link: "3_1".to_string(), 
            c_type: CType::F2,
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
            c_type: CType::F2,
            c_value: "1".to_string(),
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
        fn test_poly_h() { 
            let args = Args {
                link: "3_1".to_string(),
                c_type: CType::F2,
                c_value: "H".to_string(),
                ..Default::default()
            };
            let res = dispatch(&args);
            assert!(res.is_ok());
        }
    }
}