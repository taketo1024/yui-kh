use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use itertools::Itertools;
use yui_link::{Edge, Crossing, LinkComp};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum TngComp { 
    Arc(Edge, Edge, Edge, bool), // the middle edge keeps the min edge-id.
    Circ(Edge, bool)
}

impl TngComp { 
    pub fn short_arc(e0: Edge, e1: Edge, marked: bool) -> Self { 
        assert!(e0 != e1);
        let m = std::cmp::min(e0, e1);
        Self::arc(e0, m, e1, marked)
    }

    pub fn arc(e0: Edge, e1: Edge, e2: Edge, marked: bool) -> Self { 
        assert!(e0 != e2);
        if e0 < e2 { 
            Self::Arc(e0, e1, e2, marked)
        } else { 
            Self::Arc(e2, e1, e0, marked)
        }
    }

    pub fn circ(e: Edge, marked: bool) -> Self { 
        Self::Circ(e, marked)
    }

    pub fn from_link_comp(c: &LinkComp, base_pt: Option<Edge>) -> Self {
        let marked = base_pt.map(|e| c.edges().contains(&e)).unwrap_or(false);
        let e = c.min_edge();

        if let Some((l, r)) = c.ends() { 
            Self::arc(l, e, r, marked)
        } else { 
            Self::circ(e, marked)
        }
    }

    pub fn is_arc(&self) -> bool { 
        matches!(self, TngComp::Arc(..))
    }

    pub fn is_circle(&self) -> bool { 
        matches!(self, TngComp::Circ(..))
    }

    pub fn is_marked(&self) -> bool { 
        use TngComp::{Arc, Circ};
        match self { 
            &Arc(_, _, _, b) | &Circ(_, b) => b
        }
    }

    pub fn endpts(&self) -> Option<(Edge, Edge)> { 
        use TngComp::Arc;
        match self { 
            &Arc(e0, _, e2, _) => Some((e0, e2)),
            _ => None
        }
    }

    pub fn min_edge(&self) -> Edge { 
        use TngComp::{Arc, Circ};
        match self { 
            &Arc(_, e, _, _) | &Circ(e, _) => e
        }
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        use TngComp::Arc;
        let Arc(l0, _, l2, _) = self  else { return false };
        let Arc(r0, _, r2, _) = other else { return false };
        l0 == r0 || l0 == r2 || l2 == r0 || l2 == r2
    }

    pub fn connect(&mut self, other: Self) { 
        use TngComp::Arc;

        let Arc(l0, l1, l2, b1) = *self else { panic!() };
        let Arc(r0, r1, r2, b2) = other else { panic!() };

        let e1 = std::cmp::min(l1, r1);
        let (e0, e2) = 
        if l2 == r0 {
            (l0, r2)
        } else if l2 == r2 {
            (l0, r0)
        } else if l0 == r0 {
            (l2, r2)
        } else if l0 == r2 { 
            (l2, r0)
        } else {
            panic!()
        };
        let marked = b1 || b2;

        *self = if e0 == e2 { 
            TngComp::circ(e1, marked)
        } else { 
            TngComp::arc(e0, e1, e2, marked)
        }
    }
}

impl Display for TngComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use TngComp::{Arc, Circ};
        match self { 
            Arc(e0, _, e2, false) => write!(f, "-{e0}-{e2}-"),
            Arc(e0, _, e2,  true) => write!(f, "-{e0}-*-{e2}-"),
            Circ(e, false)        => write!(f, "○{}", yui::util::format::subscript(*e as isize)),
            Circ(e,  true)        => write!(f, "*○{}", yui::util::format::subscript(*e as isize))
        }
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Tng {
    comps: Vec<TngComp> // arc or circle
}

impl Tng { 
    pub fn new(comps: Vec<TngComp>) -> Self { 
        let mut tng = Self { comps };
        tng.comps.sort();
        tng
    }

    pub fn from_resolved(x: &Crossing, base_pt: Option<Edge>) -> Self { 
        assert!(x.is_resolved());

        let (r0, r1) = x.arcs();
        let (mut c0, c1) = (
            TngComp::from_link_comp(&r0, base_pt), 
            TngComp::from_link_comp(&r1, base_pt)
        );

        if c0.is_connectable(&c1) { 
            c0.connect(c1);
            Self::from(c0)
        } else { 
            Self::new(vec![c0, c1])
        }
    }

    pub fn empty() -> Self { 
        Self::new(vec![])
    }

    pub fn is_empty(&self) -> bool { 
        self.comps.is_empty()
    }

    pub fn is_closed(&self) -> bool { 
        self.comps.iter().all(|a| a.is_circle())
    }

    pub fn ncomps(&self) -> usize { 
        self.comps.len()
    }

    pub fn comps(&self) -> &Vec<TngComp> { 
        &self.comps
    }

    pub fn comp(&self, i: usize) -> &TngComp { 
        &self.comps[i]
    }

    pub fn sub<I>(&self, iter: I) -> Self 
    where I: IntoIterator<Item = usize> { 
        Self::new(
            iter.into_iter().map(|i| *self.comp(i)).collect()
        )
    }

    pub fn endpts(&self) -> HashSet<Edge> { 
        self.comps.iter().flat_map(|c| 
            c.endpts().map(|(e0, e1)| vec![e0, e1]).unwrap_or_default()
        ).collect()
    }

    pub fn contains(&self, c: &TngComp) -> bool { 
        self.comps.contains(c)
    }

    pub fn index_of(&self, c: &TngComp) -> Option<usize> {
        self.comps.iter().position(|c1| c1 == c)
    }

    pub fn remove_at(&mut self, i: usize) -> TngComp {
        self.comps.remove(i)
    }

    pub fn connect(&mut self, other: Self) {
        for c in other.comps.into_iter() {
            if c.is_circle() { 
                self.comps.push(c);
            } else { 
                self.append_arc(c);
            }
        }
        self.comps.sort()
    }

    pub fn append_arc(&mut self, arc: TngComp) { 
        assert!(arc.is_arc());

        let n = self.ncomps();

        // If one end of `arc` is connectable:
        if let Some(i) = self.find_connectable(&arc, n) { 
            self.comps[i].connect(arc);

            // If the other end is also connectable to a different component:
            if let Some(j) = self.find_connectable(&self.comps[i], i) { 
                let arc_j = self.comps.remove(j);
                self.comps[i].connect(arc_j);
            }
        } else { 
            self.comps.push(arc);
        }
    }

    fn find_connectable(&self, arc: &TngComp, j: usize) -> Option<usize> {
        let n = self.comps.len();

        (0..n).find(|&i| 
            i != j && self.comps[i].is_connectable(arc)
        )
    }

    pub fn connected(&self, other: &Self) -> Self { 
        let mut res = self.clone();
        res.connect(other.clone());
        res
    }

    pub fn find_loop(&self, allow_marked: bool) -> Option<usize> {
        self.comps.iter().enumerate().find(|(_, c)| 
            c.is_circle() && (allow_marked || !c.is_marked())
        ).map(|(i, _)| i)
    }

    pub fn euler_num(&self) -> isize { 
        // NOTE: χ(arc) = 1, χ(circle) = 0. 
        self.comps.iter().filter(|c| 
            c.is_arc()
        ).count() as isize
    }

    pub fn min_comp(&self) -> Option<&TngComp> { 
        self.comps.iter().min()
    }
}

impl Display for Tng {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() { 
            write!(f, "∅")
        } else { 
            let arcs = self.comps.iter().map(|a| 
                a.to_string()
            ).join(", ");
            write!(f, "{{{}}}", arcs)
        }
    }
}

impl From<TngComp> for Tng {
    fn from(c: TngComp) -> Self {
        Self::new(vec![c])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_connectable() { 
        let c0 = TngComp::short_arc(0, 1, false);
        let c1 = TngComp::short_arc(1, 2, false);
        let c2 = TngComp::short_arc(2, 3, false);
        let e  = TngComp::circ(4, false);

        assert!( c0.is_connectable(&c1));
        assert!(!c0.is_connectable(&c2));
        assert!(!c0.is_connectable(&e));

        assert!( c1.is_connectable(&c0));
        assert!( c1.is_connectable(&c2));
        assert!(!c1.is_connectable(&e));
    }

    #[test]
    fn connect_comp() { 
        let mut c0 = TngComp::short_arc(0, 1, false);
        let c1 = TngComp::short_arc(1, 2, false);
        let c2 = TngComp::short_arc(0, 2, false);

        c0.connect(c1);
        assert_eq!(c0, TngComp::short_arc(0, 2, false));

        c0.connect(c2);
        assert_eq!(c0, TngComp::circ(0, false));
    }

    #[test]
    fn append_arc() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);

        t.append_arc(TngComp::short_arc(0, 1, false)); // [0-1]
        assert_eq!(t.ncomps(), 1);

        t.append_arc(TngComp::short_arc(2, 3, false)); // [0-1] [2-3]
        assert_eq!(t.ncomps(), 2);

        t.append_arc(TngComp::short_arc(1, 2, false)); // [0-1-2-3]
        assert_eq!(t.ncomps(), 1);

        t.append_arc(TngComp::short_arc(0, 3, false)); // [-0-1-2-3-]
        assert_eq!(t.ncomps(), 1);
    }

    #[test]
    fn connect() { 
        let mut t0 = Tng::new(vec![
            TngComp::short_arc(0, 1, false),
            TngComp::short_arc(2, 3, false),
            TngComp::circ(10, false),
        ]);

        let t1 = Tng::new(vec![
            TngComp::short_arc(1, 2, false),
            TngComp::short_arc(3, 4, false),
            TngComp::circ(11, false),
        ]);

        t0.connect(t1);

        assert_eq!(t0, Tng::new(vec![
            TngComp::short_arc(0, 4, false), // [0,1,2,3,4] -> [0,4]
            TngComp::circ(10, false),
            TngComp::circ(11, false),
        ]));
    }

    #[test]
    fn deloop() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);
        assert_eq!(t.find_loop(false), None);

        t.append_arc(TngComp::short_arc(0, 1, false));
        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(false), None);

        t.append_arc(TngComp::short_arc(2, 3, false));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(false), None);

        t.append_arc(TngComp::short_arc(2, 3, false));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(false), Some(1));

        t.remove_at(1);

        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(false), None);
    }
}