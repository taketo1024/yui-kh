use std::collections::HashMap;

use yui_link::{Edge, Link};

// Involutive link
#[derive(Debug, Clone)]
pub struct InvLink { 
    link: Link,
    e_map: HashMap<Edge, Edge>
}

impl InvLink { 
    pub fn new<I>(link: Link, e_pairs: I) -> InvLink
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

        Self { link, e_map }
    }
}

#[cfg(test)]
mod tests { 
    use crate::ext::LinkExt;

    use super::*;

    #[test]
    fn link() { 
        let l = Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);
        let l = InvLink::new(l, [(1,5), (2,4)]);
        println!("{:?}", l);
    }
}