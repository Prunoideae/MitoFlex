use super::helper;

use itertools::Itertools;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::sync::mpsc::Receiver;
use std::sync::mpsc::SyncSender;

pub struct FQRead(String, String, String, String);

pub fn se_read_file(fastq1: String, sender: SyncSender<FQRead>, trim: usize) {
    let mut count = 0;
    for r in helper::read_file(fastq1.as_str())
        .lines()
        .map(|x| x.unwrap())
        .tuples::<(_, _, _, _)>()
    {
        if trim != 0 && count > trim {
            break;
        }
        let _ = sender.send(FQRead(r.0, r.1, r.2, r.3));
        count += 1;
    }
}

pub fn se_trim_read(recv: Receiver<FQRead>, start: usize, end: usize, send: SyncSender<FQRead>) {
    let len = end - start;
    for fq in recv {
        let (l1, mut l2, l3, mut l4) = (fq.0, fq.1, fq.2, fq.3);
        if start != 0 {
            l2 = l2.get(start..).unwrap().to_string();
            l4 = l4.get(start..).unwrap().to_string();
        }
        if end != 0 {
            l2 = l2.get(..len).unwrap().to_string();
            l4 = l4.get(..len).unwrap().to_string();
        }
        send.send(FQRead(l1, l2, l3, l4)).unwrap();
    }
}

pub fn se_filter_n(recv: Receiver<FQRead>, ns: usize, send: SyncSender<FQRead>) {
    for fq in recv {
        if fq.1.matches('N').count() <= ns {
            send.send(fq).unwrap();
        }
    }
}

pub fn se_filter_q(recv: Receiver<FQRead>, qv: u8, limit: f32, send: SyncSender<FQRead>) {
    for fq in recv {
        let cutoff = (fq.3.len() as f32 * limit) as usize;
        if fq.3.as_bytes().iter().filter(|&x| x <= &qv).count() < cutoff {
            send.send(fq).unwrap();
        }
    }
}

pub fn se_filter_d(recv: Receiver<FQRead>, send: SyncSender<FQRead>) {
    let mut hashset: HashSet<u64> = HashSet::new();

    for fq in recv {
        let hash = calculate_hash(&(fq.1));
        if !hashset.contains(&hash) {
            hashset.insert(hash);
            send.send(fq).unwrap();
        }
    }
}

pub fn se_write_file(recv: Receiver<FQRead>, cleanq1: String) {
    let mut clean_file = helper::write_file(cleanq1.as_str());

    for fq in recv {
        writeln!(clean_file, "{}", fq.0).unwrap();
        writeln!(clean_file, "{}", fq.1).unwrap();
        writeln!(clean_file, "{}", fq.2).unwrap();
        writeln!(clean_file, "{}", fq.3).unwrap();
    }
}

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
