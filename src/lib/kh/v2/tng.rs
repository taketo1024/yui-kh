use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use delegate::delegate;
use itertools::Itertools;
use yui_link::{Edge, Crossing, LinkComp};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TngComp { 
    path: LinkComp, 
    marked: bool
}

impl TngComp { 
    // TODO: rename to `from_path`
    pub fn from_link_comp(path: LinkComp, base_pt: Option<Edge>) -> Self {
        let marked = base_pt.map(|e| path.edges().contains(&e)).unwrap_or(false);
        Self { path, marked }
    }

    pub fn arc<I>(edges: I, marked: bool) -> Self
    where I: IntoIterator<Item = Edge> { 
        let path = LinkComp::arc(edges);
        Self { path, marked }
    }

    pub fn circ<I>(edges: I, marked: bool) -> Self
    where I: IntoIterator<Item = Edge> { 
        let path = LinkComp::circ(edges);
        Self { path, marked }
    }

    delegate! { 
        to self.path { 
            pub fn is_arc(&self) -> bool;
            pub fn is_circle(&self) -> bool;
            #[call(ends)]
            pub fn endpts(&self) -> Option<(Edge, Edge)>;
        }
    }

    pub fn is_marked(&self) -> bool { 
        self.marked
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        self.path.is_connectable(&other.path)
    }

    pub fn connect(&mut self, other: Self) { 
        self.path.connect(other.path)
    }
}

impl Display for TngComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO
        self.path.fmt(f)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Tng {
    comps: Vec<TngComp> // arc or circle
}

impl Tng { 
    pub fn new(comps: Vec<TngComp>) -> Self { 
        Self { comps }
    }

    pub fn from_resolved(x: &Crossing, base_pt: Option<Edge>) -> Self { 
        assert!(x.is_resolved());

        let (r0, r1) = x.arcs();
        let (mut c0, c1) = (
            TngComp::from_link_comp(r0, base_pt), 
            TngComp::from_link_comp(r1, base_pt)
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
        let c0 = TngComp::arc([0, 1], false);
        let c1 = TngComp::arc([1, 2], false);
        let c2 = TngComp::arc([2, 3], false);
        let e  = TngComp::circ([4], false);

        assert!( c0.is_connectable(&c1));
        assert!(!c0.is_connectable(&c2));
        assert!(!c0.is_connectable(&e));

        assert!( c1.is_connectable(&c0));
        assert!( c1.is_connectable(&c2));
        assert!(!c1.is_connectable(&e));
    }

    #[test]
    fn connect_comp() { 
        let mut c0 = TngComp::arc([0, 1], false);
        let c1 = TngComp::arc([1, 2], false);
        let c2 = TngComp::arc([0, 2], false);

        c0.connect(c1);
        assert_eq!(c0, TngComp::arc([0, 1, 2], false));

        c0.connect(c2);
        assert_eq!(c0, TngComp::circ([0, 1, 2], false));
    }

    #[test]
    fn append_arc() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);

        t.append_arc(TngComp::arc([0, 1], false)); // [0-1]
        assert_eq!(t.ncomps(), 1);
        assert!(t.comp(0).is_arc());

        t.append_arc(TngComp::arc([2, 3], false)); // [0-1] [2-3]
        assert_eq!(t.ncomps(), 2);
        assert!(t.comp(0).is_arc());
        assert!(t.comp(1).is_arc());

        t.append_arc(TngComp::arc([1, 2], false)); // [0-1-2-3]
        assert_eq!(t.ncomps(), 1);
        assert!(t.comp(0).is_arc());

        t.append_arc(TngComp::arc([0, 3], false)); // [-0-1-2-3-]
        assert_eq!(t.ncomps(), 1);
        assert!(t.comp(0).is_circle());
    }

    #[test]
    fn connect() { 
        let mut t0 = Tng::new(vec![
            TngComp::arc([0, 1], false),
            TngComp::arc([2, 3], false),
            TngComp::circ([10], false),
        ]);

        let t1 = Tng::new(vec![
            TngComp::arc([1, 2], false),
            TngComp::arc([3, 4], false),
            TngComp::circ([11], false),
        ]);

        t0.connect(t1);

        assert_eq!(t0, Tng::new(vec![
            TngComp::arc([0, 1, 2, 3, 4], false),
            TngComp::circ([10], false),
            TngComp::circ([11], false),
        ]));
    }

    #[test]
    fn find_loop() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);
        assert_eq!(t.find_loop(false), None);

        t.append_arc(TngComp::arc([0, 1], false));
        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(false), None);

        t.append_arc(TngComp::arc([2, 3], false));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(false), None);

        t.append_arc(TngComp::arc([2, 3], false));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(false), Some(1));

        t.remove_at(1);

        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(false), None);
    }
}