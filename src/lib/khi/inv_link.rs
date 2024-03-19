use std::collections::HashMap;

use itertools::Itertools;
use yui_link::{Edge, Link};

// Involutive link
#[derive(Debug, Clone)]
pub struct InvLink { 
    link: Link,
    base_pt: Option<Edge>,
    e_map: HashMap<Edge, Edge>,
    x_map: HashMap<usize, usize>
}

impl InvLink { 
    pub fn new<I>(link: Link, e_pairs: I, base_pt: Option<Edge>) -> InvLink
    where I: IntoIterator<Item = (Edge, Edge)> { 
        let mut edges = link.edges();
        let mut e_map = HashMap::new();

        for (e1, e2) in e_pairs { 
            assert!(edges.contains(&e1));
            assert!(edges.contains(&e2));

            e_map.insert(e1, e2);
            e_map.insert(e2, e1);
            edges.remove(&e1);
            edges.remove(&e2);
        }

        for e in edges { 
            e_map.insert(e, e);
        }

        let mut x_map = HashMap::new();

        // TODO? check resolution
        for (i, x) in link.data().iter().enumerate() { 
            let edges = x.edges().map(|e| e_map.get(&e).unwrap());
            let find = link.data().iter().find_position(|y|
                edges.iter().all(|e| y.edges().contains(e))
            );

            assert!(find.is_some(), "no match for x: {x} -> {edges:?}");

            let j = find.unwrap().0;
            x_map.insert(i, j);
            x_map.insert(j, i);
        }

        Self { link, base_pt, e_map, x_map }
    }

    pub fn link(&self) -> &Link { 
        &self.link
    }

    pub fn base_pt(&self) -> Option<Edge> {
        self.base_pt
    }

    pub fn inv_e(&self, e: Edge) -> Edge { 
        self.e_map.get(&e).cloned().unwrap()
    }

    pub fn inv_x(&self, i: usize) -> usize { 
        self.x_map.get(&i).cloned().unwrap()
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn inv_e() { 
        let l = Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);
        let l = InvLink::new(l, [(1,5), (2,4)], None);

        assert_eq!(l.inv_e(1), 5);
        assert_eq!(l.inv_e(2), 4);
        assert_eq!(l.inv_e(3), 3);
        assert_eq!(l.inv_e(4), 2);
        assert_eq!(l.inv_e(5), 1);
        assert_eq!(l.inv_e(6), 6);
    }
    
    #[test]
    fn inv_x() { 
        let l = Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);
        let l = InvLink::new(l, [(1,5), (2,4)], None);

        assert_eq!(l.inv_x(0), 0);
        assert_eq!(l.inv_x(1), 2);
        assert_eq!(l.inv_x(2), 1);
    }
}