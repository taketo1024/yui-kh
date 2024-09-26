use std::collections::{HashMap, HashSet};
use std::fmt::Display;

use log::info; 
use itertools::Itertools;
use num_traits::Zero;
use cartesian::cartesian;
use yui::{Ring, RingOps, Sign};
use yui_homology::{XChainComplex, XModStr, Grid1};
use yui_link::{Crossing, Edge, State};
use yui::bitseq::Bit;

use crate::kh::{KhAlgGen, KhGen, KhLabel};
use super::cob::{Cob, Dot, Bottom, CobComp, LcCob, LcCobTrait};
use super::tng::{Tng, TngComp};

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct TngKey { 
    pub(crate) state: State,
    pub(crate) label: KhLabel
}

impl TngKey { 
    pub(crate) fn init() -> Self { 
        Self { state: State::empty(), label: KhLabel::empty() }
    }

    fn as_gen(&self, deg_shift: (isize, isize)) -> KhGen { 
        KhGen::new(self.state, self.label, deg_shift)
    }

    fn append(&mut self, other: TngKey) { 
        self.state.append(other.state);
        self.label.append(other.label);
    }

    fn appended(&self, other: &TngKey) -> TngKey { 
        let mut res = self.clone();
        res.append(other.clone());
        res
    }

    fn add_label(&self, l: KhAlgGen) -> TngKey {
        let mut res = self.clone();
        res.label.push(l);
        res
    }
}

impl From<&KhGen> for TngKey {
    fn from(x: &KhGen) -> Self {
        TngKey { state: x.state, label: x.label }
    }
}

impl Display for TngKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.as_gen((0, 0)).fmt(f)
    }
}

#[derive(Clone, Debug)]
pub struct TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    key: TngKey,
    tng: Tng,
    in_edges: HashSet<TngKey>,
    out_edges: HashMap<TngKey, LcCob<R>>
}

impl<R> TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn init() -> Self { 
        let key = TngKey::init();
        let tng = Tng::empty();
        let in_edges = HashSet::new();
        let out_edges = HashMap::new();
        Self { key, tng, in_edges, out_edges }
    }

    pub fn tng(&self) -> &Tng { 
        &self.tng
    }

    pub(crate) fn out_edges(&self) -> &HashMap<TngKey, LcCob<R>> {
        &self.out_edges
    }

    fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        let key = self.key;
        let tng = self.tng.convert_edges(&f);
        let in_edges = self.in_edges.clone();
        let out_edges = self.out_edges.iter().map(|(k, cob)|
            (k.clone(), cob.convert_edges(&f))
        ).collect();
        TngVertex { key, tng, in_edges, out_edges }
    }
}

impl<R> Display for TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}; {})", self.key, self.tng)
    }
}

#[derive(Debug, Default)]
pub struct TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    h: R,
    t: R,
    deg_shift: (isize, isize),
    base_pt: Option<Edge>,
    vertices: HashMap<TngKey, TngVertex<R>>,
    crossings: Vec<Crossing>,
}

impl<R> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>, vertices: HashMap<TngKey, TngVertex<R>>, crossings: Vec<Crossing>) -> Self { 
        let (h, t) = (h.clone(), t.clone());
        TngComplex{ h, t, deg_shift, base_pt, vertices, crossings }
    }

    pub fn init(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let mut vertices = HashMap::new();
        let k0 = TngKey::init();
        let v0 = TngVertex::init();
        vertices.insert(k0, v0);

        TngComplex::new(h, t, deg_shift, base_pt, vertices, vec![])
    }

    pub fn from_crossing(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>, x: &Crossing) -> Self { 
        let mut _self = Self::new(h, t, deg_shift, base_pt, HashMap::new(), vec![]);

        if x.is_resolved() { 
            let mut v = TngVertex::init();
            v.tng = Tng::from_resolved(x);
            _self.vertices.insert(v.key, v);
        } else { 
            let mut v0 = TngVertex::init();
            v0.key.state.push_0();
            v0.tng = Tng::from_resolved(&x.resolved(Bit::Bit0));
            let k0 = v0.key;

            let mut v1 = TngVertex::init();
            v1.key.state.push_1();
            v1.tng = Tng::from_resolved(&x.resolved(Bit::Bit1));
            let k1 = v1.key;

            _self.vertices.insert(k0, v0);
            _self.vertices.insert(k1, v1);

            let sdl = LcCob::from(Cob::from(CobComp::sdl_from(&x)));
            _self.add_edge(&k0, &k1, sdl);

            _self.crossings.push(x.clone());
        }

        _self
    }

    pub fn ht(&self) -> (&R, &R) { 
        (&self.h, &self.t)
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn base_pt(&self) -> Option<Edge> { 
        self.base_pt
    }

    pub fn dim(&self) -> usize { 
        self.crossings.iter().filter(|x| !x.is_resolved()).count()
    }

    pub fn rank(&self, i: isize) -> usize { 
        let i0 = self.deg_shift.0;
        self.vertices.keys().filter(|k| {
            let w = k.state.weight() as isize;
            w == i - i0
        }).count()
    }

    pub fn crossing(&self, i: usize) -> &Crossing {
        &self.crossings[i]
    }

    pub fn find_crossing(&self, x: &Crossing) -> Option<usize> { 
        self.crossings.iter().position(|y| x == y)
    }

    pub fn vertex(&self, v: &TngKey) -> &TngVertex<R> { 
        &self.vertices[v]
    }

    pub fn nverts(&self) -> usize { 
        self.vertices.len()
    }

    pub fn iter_verts(&self) -> impl Iterator<Item = (&TngKey, &TngVertex<R>)> {
        self.vertices.iter().sorted_by(|(k0, _), (k1, _)| k0.cmp(k1))
    }

    fn keys_into(&self, k: &TngKey) -> impl Iterator<Item = &TngKey> { 
        self.vertex(k).in_edges.iter()
    }

    fn keys_out_from(&self, k: &TngKey) -> impl Iterator<Item = &TngKey> { 
        self.vertex(k).out_edges.keys()
    }

    fn remove_vertex(&mut self, k: &TngKey) { 
        let in_edges = self.keys_into(k).cloned().collect_vec();
        let out_edges = self.keys_out_from(k).cloned().collect_vec();

        self.vertices.remove(k);

        for j in in_edges { 
            self.vertices.get_mut(&j).unwrap().out_edges.remove(k);
        }
        
        for l in out_edges { 
            self.vertices.get_mut(&l).unwrap().in_edges.remove(k);
        }
    }
    
    fn rename_vertex_key(&mut self, k_old: &TngKey, k_new: TngKey) { 
        assert_ne!(k_old, &k_new);

        let in_edges = self.keys_into(k_old).cloned().collect_vec(); 
        let in_removed = in_edges.into_iter().map(|j| {
            let f = self.remove_edge(&j, k_old);
            (j, f)
        }).collect_vec();

        let out_edges = self.keys_out_from(k_old).cloned().collect_vec();
        let out_removed = out_edges.into_iter().map(|l| {
            let f = self.remove_edge(k_old, &l);
            (l, f)
        }).collect_vec();

        let mut v = self.vertices.remove(k_old).unwrap();
        v.key = k_new;
        self.vertices.insert(k_new, v);

        for (j, f) in in_removed { 
            self.add_edge(&j, &k_new, f);
        }

        for (l, f) in out_removed { 
            self.add_edge(&k_new, &l, f);
        }
    }

    fn duplicate_vertex(&mut self, k: &TngKey, k_new: TngKey) { 
        assert_ne!(k, &k_new);

        let in_edges = self.vertex(k).in_edges.clone(); 
        let out_edges = self.vertex(k).out_edges.keys().cloned().collect_vec();

        let mut v_new = self.vertex(k).clone();
        v_new.key = k_new;
        v_new.in_edges.clear();
        v_new.out_edges.clear();

        self.vertices.insert(k_new, v_new);

        for j in in_edges { 
            let f = self.edge(&j, k).clone();
            self.add_edge(&j, &k_new, f);
        }

        for l in out_edges { 
            let f = self.edge(k, &l).clone();
            self.add_edge(&k_new, &l, f);
        }
    }
    
    pub fn edge(&self, k: &TngKey, l: &TngKey) -> &LcCob<R> {
        &self.vertices[k].out_edges[l]
    }

    fn has_edge(&self, k: &TngKey, l: &TngKey) -> bool { 
        debug_assert_eq!(
            self.vertices[k].out_edges.contains_key(l), 
            self.vertices[l].in_edges.contains(k)
        );
        self.vertices[k].out_edges.contains_key(l)
    }

    fn add_edge(&mut self, k: &TngKey, l: &TngKey, f: LcCob<R>) { 
        assert!(!self.has_edge(k, l));
        assert!(!f.is_zero());

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.insert(*l, f);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.insert(*k);
    }

    fn remove_edge(&mut self, k: &TngKey, l: &TngKey) -> LcCob<R> { 
        assert!(self.has_edge(k, l));
        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.remove(k);

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.remove(l).unwrap()
    }

    fn modify_edge<F>(&mut self, k: &TngKey, l: &TngKey, map: F)
    where F: Fn(LcCob<R>) -> LcCob<R> {
        assert!(self.has_edge(k, l));

        let f = self.remove_edge(k, l);
        let map_f = map(f);

        if !map_f.is_zero() { 
            self.add_edge(k, l, map_f);
        }
    }
    
    // TODO rename to add_crossing
    pub fn append(&mut self, x: &Crossing) {
        info!("({}) append: {x}", self.nverts());

        let (h, t) = self.ht();
        let c = Self::from_crossing(h, t, (0, 0), None, x);
        self.connect(c);
    }

    // See [Bar-Natan '05] Section 5.
    // https://arxiv.org/abs/math/0410495
    pub fn connect(&mut self, mut other: TngComplex<R>) { 
        assert_eq!(self.ht(), other.ht());
        assert!(self.base_pt.is_none() || other.base_pt.is_none() || self.base_pt == other.base_pt);

        // create vertices
        let vertices = std::mem::take(&mut self.vertices);
        self.vertices = cartesian!(vertices.iter(), other.vertices.iter()).map(|((k, v), (l, w))| {  
            let kl = k.appended(l);
            let mut vw = TngVertex::init();
            vw.key = kl;
            vw.tng = v.tng.connected(&w.tng); // D(v, w)
            (kl, vw)
        }).collect();

        // create edges
        cartesian!(vertices.iter(), other.vertices.iter()).for_each(|((k0, v0), (l0, w0))| {
            let i0 = (k0.state.weight() as isize) - self.deg_shift.0;
            let k0_l0 = k0.appended(l0);
            
            for (k1, f) in v0.out_edges.iter() { 
                let k1_l0 = k1.appended(l0);
                let f_id = f.connected(&Cob::id(w0.tng())); // D(f, 1) 
                
                if !f_id.is_zero() { 
                    self.add_edge(&k0_l0, &k1_l0, f_id);
                }                
            }

            for (l1, f) in w0.out_edges.iter() { 
                let k0_l1 = k0.appended(l1);
                let e = R::from_sign(Sign::from_parity(i0 as i64));
                let id_f = f.connected(&Cob::id(v0.tng())) * e; // (-1)^{deg(k0)} D(1, f) 
                if !id_f.is_zero() { 
                    self.add_edge(&k0_l0, &k0_l1, id_f);
                }
            }
        });

        self.base_pt = self.base_pt.or(other.base_pt);
        self.deg_shift.0 += other.deg_shift.0;
        self.deg_shift.1 += other.deg_shift.1;
        self.crossings.append(&mut other.crossings);

        // self.validate();
    }

    pub fn find_comp<F>(&self, pred: F) -> Option<(TngKey, usize)>
    where F: Fn(&TngComp) -> bool {
        for (k, v) in self.iter_verts() { 
            if let Some(r) = v.tng.find_comp(&pred) { 
                return Some((*k, r))
            }
        }
        None
    }

    pub fn find_loop(&self, allow_marked: bool) -> Option<(TngKey, usize)> { 
        self.find_comp(|c|
            c.is_circle() && (allow_marked || self.base_pt.map(|e| !c.contains(e)).unwrap_or(true))
        )
    }

    pub fn find_merge_or_split_loop(&self, allow_marked: bool) -> Option<(TngKey, usize)> { 
        let base_pt = self.base_pt;
        let ok = |c: &TngComp| {
            c.is_circle() && (allow_marked || base_pt.map(|e| !c.contains(e)).unwrap_or(true))
        };
        
        self.vertices.iter().find_map(|(k, v)|
            v.out_edges.iter().find_map(|(l, f)|
                f.gens().find_map(|cob| 
                    cob.comps().find_map(|cob| 
                        if cob.is_merge() { 
                            let Some(r_loc) = cob.src().find_comp(&ok) else { 
                                return None
                            };
                            let circ = cob.src().comp(r_loc);
                            let r = v.tng.find_comp(|c| c == circ).unwrap();
                            Some((*k, r))
                        } else if cob.is_split() {
                            let Some(r_loc) = cob.tgt().find_comp(&ok) else { 
                                return None
                            };
                            let circ = cob.tgt().comp(r_loc);
                            let r = self.vertex(l).tng.find_comp(|c| c == circ).unwrap();
                            Some((*l, r))
                        } else { 
                            None
                        }
                    )
                )
            )
        )
    }

    pub fn deloop(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> { 
        info!("({}) deloop {} at {r}", self.nverts(), &self.vertices[k]);

        let c = self.vertex(k).tng.comp(r);
        assert!(c.is_circle());
        
        let marked = self.base_pt.map(|e| c.contains(e)).unwrap_or(false);

        #[allow(non_snake_case)]
        let updated_keys = if marked { 
            let k_X = k.add_label(KhAlgGen::X);

            self.rename_vertex_key(k, k_X);
            self.deloop_with(&k_X, r, Dot::X, Dot::None);

            vec![k_X]
        } else { 
            let k_X = k.add_label(KhAlgGen::X);
            let k_1 = k.add_label(KhAlgGen::I);

            self.rename_vertex_key(k, k_X);
            self.duplicate_vertex(&k_X, k_1);

            self.deloop_with(&k_X, r, Dot::X, Dot::None);
            self.deloop_with(&k_1, r, Dot::None, Dot::Y);

            vec![k_X, k_1]
        };

        updated_keys
    }

    fn deloop_with(&mut self, k: &TngKey, r: usize, birth_dot: Dot, death_dot: Dot) { 
        // remove circle
        let circ = self.vertices.get_mut(k).unwrap().tng.remove_at(r);

        let v_in = self.keys_into(k).cloned().collect_vec();
        let v_out = self.keys_out_from(k).cloned().collect_vec();

        // cap incoming cobs
        let (h, t) = (self.h.clone(), self.t.clone());
        for j in v_in.iter() { 
            self.modify_edge(j, k, |f|
                f.cap_off(Bottom::Tgt, &circ, death_dot).part_eval(&h, &t)
            );
        }
        
        // cup outgoing cobs
        for l in v_out.iter() { 
            self.modify_edge(k, l, |f|
                f.cap_off(Bottom::Src, &circ, birth_dot).part_eval(&h, &t)
            );
        }
    }

    pub fn find_edge<'a, 'b, F>(&'a self, k: &'a TngKey, pred: F) -> Option<(TngKey, TngKey)>
    where 'a: 'b, F: Fn(&'b TngKey, &'b TngKey) -> bool { 
        let mut cand = None;
        let mut cand_s = usize::MAX;

        // collect candidate edges into and out from v. 

        let v = self.vertex(k);
        let edges = Iterator::chain(
            v.in_edges.iter().map(|i| (i, k)),
            v.out_edges.keys().map(|i| (k, i))
        );

        let score = |i, j| {
            let ni = self.vertex(i).out_edges.len(); // nnz in column i
            let nj = self.vertex(j).in_edges.len();  // nnz in row j
            (ni - 1) * (nj - 1)
        };

        for (i, j) in edges { 
            if pred(i, j) { 
                let s = score(i, j);

                if s == 0 {
                    return Some((*i, *j));
                } else if s < cand_s { 
                    cand = Some((i, j));
                    cand_s = s;
                }
            }
        }

        if let Some((i, j)) = cand { 
            Some((*i, *j))
        } else { 
            None
        }
    }

    pub fn find_inv_edge(&self, k: &TngKey) -> Option<(TngKey, TngKey)> { 
        let pred = |i, j| self.edge(i, j).is_invertible();
        self.find_edge(k, pred)
    }

    //  Gaussian Elimination
    //
    //       a
    //  v0 - - -> v1         .             .
    //     \   / b
    //       /         ==>  
    //     /   \ c              d - ca⁻¹b
    //  w0 -----> w1         w0 ---------> w1
    //       d                

    pub fn eliminate(&mut self, k0: &TngKey, k1: &TngKey) {
        let a = self.edge(k0, k1);
        let Some(ainv) = a.inv() else { 
            panic!("{a} is not invertible.")
        };

        info!("({}) eliminate {}: {} -> {}", self.nverts(), a, self.vertex(k0), self.vertex(k1));
        
        let v0_out = self.keys_out_from(k0).cloned().collect_vec();
        let v1_in  = self.keys_into(k1).cloned().collect_vec();
        
        for (l0, l1) in cartesian!(v1_in.iter(), v0_out.iter()) {
            if l0 == k0 || l1 == k1 { 
                continue
            }

            let (h, t) = (&self.h, &self.t);
            
            let b = self.edge(l0, k1);
            let c = self.edge(k0, l1);
            let cab = (c * &ainv * b).part_eval(h, t);

            if cab.is_zero() { 
                continue
            } 

            let s = if self.has_edge(l0, l1) { 
                let d = self.remove_edge(l0, l1);
                d - cab
            } else { 
                -cab
            };

            if !s.is_zero() { 
                self.add_edge(l0, l1, s);
            }
        }

        self.remove_vertex(k0);
        self.remove_vertex(k1);

        // self.validate();
    }

    pub fn into_complex(self) -> XChainComplex<KhGen, R> {
        debug_assert!(self.is_finalizable());

        let n = self.dim();
        let i0 = self.deg_shift.0;
        let i1 = i0 + (n as isize);

        let summands = Grid1::generate(i0..=i1, |i| { 
            let w = (i - i0) as usize;
            let gens = self.vertices.keys().filter_map(|k| 
                (k.state.weight() == w).then(|| 
                    k.as_gen(self.deg_shift)
                )
            ).sorted();
            XModStr::free(gens)
        });

        XChainComplex::new(summands, 1, move |_, z| { 
            let (h, t) = self.ht();
            z.apply(|x| {
                let k = TngKey::from(x);
                let v = self.vertex(&k);
                v.out_edges.iter().map(|(l, f)|
                    (l.as_gen(x.deg_shift), f.eval(h, t))
                ).collect()
            })
        })
    }

    fn is_finalizable(&self) -> bool { 
        self.vertices.iter().all(|(_, v)|
            v.tng.is_empty()
        )
    }

    pub fn desc_d(&self) -> String { 
        let mut str = "".to_string();
        for k0 in self.vertices.keys().sorted() { 
            let v = &self.vertices[&k0];
            str += &format!("{k0}: {}", v.tng);

            for k1 in v.out_edges.keys().sorted() { 
                let f = &v.out_edges[&k1];
                str += &format!("\n  -> {k1}: {f}");
            }
            str += "\n";
        }
        str
    }

    pub fn print_d(&self) { 
        println!("{}", self.desc_d());
    }

    pub fn validate(&self) {
        for (k, v) in self.vertices.iter() { 
            // validate in_edges 
            for j in self.keys_into(k) {
                assert!(
                    self.vertices.contains_key(j),
                    "no vertex for in-edge {j} -> {k}"
                );

                let u = self.vertex(j);
                
                assert!(
                    u.out_edges.contains_key(k),
                    "no out-edge {j} -> {k}"
                );
            }
            
            // validate out_edges 
            for l in self.keys_out_from(k) {
                assert!(
                    self.vertices.contains_key(l),
                    "no vertex for out-edge {k} -> {l}"
                );

                let w = self.vertex(l);

                assert!(
                    w.in_edges.contains(k),
                    "no in-edge {k} -> {l}"
                );
            }

            // validate cobordism
            for l in self.keys_out_from(k) {
                let w = self.vertex(l);
                let f = self.edge(k, l);

                assert!(!f.is_zero());

                f.iter().for_each(|(cob, _)| { 
                    assert_eq!(&cob.src(), v.tng(), "invalid source: {} for {cob}", v.tng());
                    assert_eq!(&cob.tgt(), w.tng(), "invalid target: {} for {cob}", w.tng());
                })
            }
        }
    }

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        let base_pt = self.base_pt.map(|e| f(e));
        let crossings = self.crossings.iter().map(|x| x.convert_edges(&f)).collect();

        let vertices = self.iter_verts().map(|(k1, v1)| {
            let k2 = k1.clone();
            let v2 = v1.convert_edges(&f);
            (k2, v2)
        }).collect();

        TngComplex::new(&self.h, &self.t, self.deg_shift, base_pt, vertices, crossings)
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use crate::kh::KhLabel;

    #[test]
    fn empty() { 
        let c = TngComplex::init(&0, &0, (0, 0), None);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn single_x() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x = Crossing::from_pd_code([0,1,2,3]);
        c.append(&x);

        assert_eq!(c.dim(), 1);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
    }

    #[test]
    fn single_x_resolved() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x = Crossing::from_pd_code([0,1,2,3]).resolved(Bit::Bit0);
        c.append(&x);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn two_x_disj() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([0,1,2,3]);
        let x1 = Crossing::from_pd_code([4,5,6,7]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn two_x() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([0,4,1,5]);
        let x1 = Crossing::from_pd_code([3,1,4,2]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn deloop() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([0, 1, 1, 0]).resolved(Bit::Bit0); // unknot
        c.append(&x0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        let e = c.find_loop(false);
        assert!(e.is_some());
        
        let Some((k, r)) = e else { panic!() };
        let keys = c.deloop(&k, r);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 2);
        assert_eq!(keys.len(), 2);
    }

    #[test]
    fn deloop_tangle() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([4,2,5,1]);
        let x1 = Crossing::from_pd_code([3,6,4,1]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        let e = c.find_loop(false);
        assert!(e.is_some());

        let Some((k, r)) = e else { panic!() };

        assert_eq!(k, TngKey {
            state: State::from([1,0]), 
            label: KhLabel::from_iter([])
        });
        assert_eq!(r, 2);

        c.deloop(&k, r);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 3); // delooped here
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn deloop_marked() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), Some(0)); // base point = 0
        let x0 = Crossing::from_pd_code([0, 1, 1, 0]).resolved(Bit::Bit0); // unknot
        c.append(&x0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        let e = c.find_loop(false);
        assert!(e.is_none());
        
        let e = c.find_loop(true);
        assert!(e.is_some());
        
        let Some((k, r)) = e else { panic!() };
        let keys = c.deloop(&k, r);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1); // not 2
        assert_eq!(keys.len(), 1); // not 2
    }

    #[test]
    fn connect() {
        let mut c0 = TngComplex::init(&0, &0, (0, 0), None);
        let mut c1 = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([4,2,5,1]);
        let x1 = Crossing::from_pd_code([3,6,4,1]);

        c0.append(&x0);
        c1.append(&x1);

        c0.connect(c1);

        assert_eq!(c0.dim(), 2);
        assert_eq!(c0.rank(0), 1);
        assert_eq!(c0.rank(1), 2);
        assert_eq!(c0.rank(2), 1);

        c0.validate();
    }

    #[test]
    fn connect_trefoil() {
        let mut c0 = TngComplex::init(&0, &0, (0, 0), None);
        let mut c1 = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([1,4,2,5]);
        let x1 = Crossing::from_pd_code([3,6,4,1]);
        let x2 = Crossing::from_pd_code([5,2,6,3]);

        c0.append(&x0);
        c0.append(&x1);
        c1.append(&x2);

        c0.connect(c1);

        assert_eq!(c0.dim(), 3);
        assert_eq!(c0.rank(0), 1);
        assert_eq!(c0.rank(1), 3);
        assert_eq!(c0.rank(2), 3);
        assert_eq!(c0.rank(3), 1);

        c0.validate();
    }
}
