use super::helper;

use itertools::Itertools;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::sync::mpsc::Receiver;
use std::sync::mpsc::SyncSender;

pub struct FQRead2(
    (String, String, String, String),
    (String, String, String, String),
);

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
