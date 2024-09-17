use std::marker::PhantomData;
use std::str::FromStr;
use yui::{EucRing, EucRingOps};
use yui_homology::{DisplayForGrid, DisplaySeq, DisplayTable, GridDeg, GridTrait, RModStr, XHomologyBase, XModStr};
use yui_kh::kh::{KhChain, KhGen};
use yui_kh::khi::{KhIChain, KhIChainExt, KhIComplex, KhIGen, KhIHomology};
use yui_link::Link;
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
        ensure!(t.is_zero(), "");
    
        let l = load_sinv_knot(&self.args.link, self.args.mirror)?;
        let ckhi = KhIComplex::<R>::new(&l, &h, self.args.reduced);

        let poly = ["H", "0,T"].contains(&self.args.c_value.as_str());
        let bigraded = h.is_zero() && t.is_zero();

        if bigraded { 
            let with_trans = self.args.show_gens || self.args.show_alpha;

            let ckhi = ckhi.into_bigraded();
            let khi = ckhi.homology(with_trans);
    
            self.out(&khi.display_table("i", "j"));

            if self.args.show_gens { 
                self.show_gens(&khi);
            }
        } else if poly { 
            let with_trans = true;

            let khi = ckhi.homology(with_trans);
            let gens = khi.gen_table();
    
            self.out(&gens.display_table("i", "j"));
    
            if self.args.show_gens { 
                self.show_gens(khi.inner());
            }

            let zs = ckhi.canon_cycles();

            if self.args.show_alpha { 
                self.show_alpha(&khi, zs);
            }
        } else { 
            let with_trans = self.args.show_gens || self.args.show_alpha;

            let khi = ckhi.homology(with_trans);
            self.out(&khi.display_seq("i"));

            if self.args.show_gens { 
                self.show_gens(khi.inner());
            }

            let zs = ckhi.canon_cycles();

            if self.args.show_alpha { 
                self.show_alpha(&khi, zs);
            }
        }

        Ok(self.flush())
    }

    fn show_gens<I>(&mut self, kh: &XHomologyBase<I, KhIGen, R>)
    where I: GridDeg { 
        for i in kh.support() {
            let h = &kh[i];
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

    fn show_ssi(&mut self, l: &Link, c: &R, kh0: &XModStr<KhGen, R>, zs: &[KhChain<R>]) { 
        assert!(!c.is_unit() && !c.is_unit());

        use yui_kh::misc::div_vec;
        let Some(a) = zs.first() else { return };

        let w = l.writhe();
        let r = l.seifert_circles().len() as i32;
        let v = kh0.vectorize(a).subvec(0..kh0.rank());
        let d = div_vec(&v, c).unwrap();
        let s = 2 * d + w - r + 1;

        self.out(&format!("ss = {s} (d = {d}, w = {w}, r = {r})"));
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