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

pub fn pe_read_file(fastq1: String, fastq2: String, sender: SyncSender<FQRead2>, trim: usize) {
    let mut count = 0;
    for (r1, r2) in helper::read_file(fastq1.as_str())
        .lines()
        .map(|x| x.unwrap())
        .tuples::<(_, _, _, _)>()
        .zip(
            helper::read_file(fastq2.as_str())
                .lines()
                .map(|x| x.unwrap())
                .tuples::<(_, _, _, _)>(),
        )
    {
        if trim != 0 && count > trim {
            break;
        }
        sender
            .send(FQRead2((r1.0, r1.1, r1.2, r1.3), (r2.0, r2.1, r2.2, r2.3)))
            .unwrap();
        count += 1;
    }
}

pub fn pe_trim_read(
    recv: Receiver<FQRead2>,
    start: usize,
    end: usize,
    sender: SyncSender<FQRead2>,
) {
    let len = end - start;
    for FQRead2(fq1, fq2) in recv {
        let (l11, mut l21, l31, mut l41) = (fq1.0, fq1.1, fq1.2, fq1.3);
        let (l12, mut l22, l32, mut l42) = (fq2.0, fq2.1, fq2.2, fq2.3);
        if start != 0 {
            l21 = l21.get(start..).unwrap().to_string();
            l41 = l41.get(start..).unwrap().to_string();

            l22 = l22.get(start..).unwrap().to_string();
            l42 = l42.get(start..).unwrap().to_string();
        }

        if end != 0 {
            l21 = l21.get(..len).unwrap().to_string();
            l41 = l41.get(..len).unwrap().to_string();

            l22 = l22.get(..len).unwrap().to_string();
            l42 = l42.get(..len).unwrap().to_string();
        }

        sender
            .send(FQRead2((l11, l21, l31, l41), (l12, l22, l32, l42)))
            .unwrap();
    }
}

pub fn pe_filter_n(recv: Receiver<FQRead2>, ns: usize, sender: SyncSender<FQRead2>) {
    for fqr2 in recv {
        let fq1 = &fqr2.0;
        let fq2 = &fqr2.1;

        if fq1.1.matches('N').count() <= ns && fq2.1.matches('N').count() <= ns {
            sender.send(fqr2).unwrap();
        }
    }
}

pub fn pe_filter_q(recv: Receiver<FQRead2>, qv: u8, limit: f32, sender: SyncSender<FQRead2>) {
    for fqr2 in recv {
        let fq1 = &fqr2.0;
        let fq2 = &fqr2.1;

        let len1 = (fq1.3.len() as f32 * limit) as usize;
        let len2 = (fq2.3.len() as f32 * limit) as usize;

        if fq1.3.as_bytes().iter().filter(|&x| x <= &qv).count() <= len1
            && fq2.3.as_bytes().iter().filter(|&x| x <= &qv).count() <= len2
        {
            sender.send(fqr2).unwrap();
        }
    }
}

pub fn pe_filter_d(recv: Receiver<FQRead2>, sender: SyncSender<FQRead2>) {
    let mut dup1: HashSet<u64> = HashSet::new();
    let mut dup2: HashSet<u64> = HashSet::new();

    for fqr2 in recv {
        let fq1 = &(fqr2.0);
        let fq2 = &(fqr2.1);

        let h1 = calculate_hash(&(fq1.1));
        let h2 = calculate_hash(&(fq2.1));

        if !dup1.contains(&h1) || !dup2.contains(&h2) {
            sender.send(fqr2).unwrap();
            dup1.insert(h1);
            dup2.insert(h2);
        }
    }
}

pub fn pe_write_file(recv: Receiver<FQRead2>, cleanq1: String, cleanq2: String) {
    let mut cleanf1 = helper::write_file(cleanq1.as_str());
    let mut cleanf2 = helper::write_file(cleanq2.as_str());

    for FQRead2(r1, r2) in recv {
        writeln!(cleanf1, "{}", r1.0).unwrap();
        writeln!(cleanf1, "{}", r1.1).unwrap();
        writeln!(cleanf1, "{}", r1.2).unwrap();
        writeln!(cleanf1, "{}", r1.3).unwrap();

        writeln!(cleanf2, "{}", r2.0).unwrap();
        writeln!(cleanf2, "{}", r2.1).unwrap();
        writeln!(cleanf2, "{}", r2.2).unwrap();
        writeln!(cleanf2, "{}", r2.3).unwrap();
    }
}

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
