use super::helper;

use itertools::Itertools;
use rayon::ThreadPool;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::sync::mpsc::Receiver;
use std::sync::mpsc::SyncSender;
use std::sync::Arc;

pub struct FQRead(String, String, String, String);

pub fn se_reader_worker(
    fastq: String,
    trim: usize,
    cpupool: ThreadPool,
    sender: Arc<SyncSender<Vec<FQRead>>>,
    ns: usize,
    qv: u8,
    limit: f32,
    start: usize,
    end: usize,
) {
    let mut count = 0;
    let mut chunked = Vec::with_capacity(10000);

    for fq in helper::read_file(fastq.as_str())
        .lines()
        .map(|x| x.unwrap())
        .tuples::<(_, _, _, _)>()
    {
        chunked.push(FQRead(fq.0, fq.1, fq.2, fq.3));
        if chunked.len() >= 10000 {
            let sender_child = sender.clone();
            cpupool.spawn(move || {
                se_calculate_unit(chunked, &sender_child, ns, qv, limit, start, end)
            });
            chunked = Vec::with_capacity(10000);
        }
        count += 1;
        if count > trim {
            break;
        }
    }
}

fn se_calculate_unit(
    chunked_reads: Vec<FQRead>,
    sender: &SyncSender<Vec<FQRead>>,
    ns: usize,
    qv: u8,
    limit: f32,
    start: usize,
    end: usize,
) {
    let mut chunked_reads = chunked_reads;

    if start != 0 || end != 0 {
        let len = end - start;
        chunked_reads = chunked_reads
            .into_iter()
            .map(|x| {
                let mut x2 = x.1;
                let mut x4 = x.3;
                if start != 0 {
                    x2 = x2.get(start..).unwrap().to_string();
                    x4 = x4.get(start..).unwrap().to_string();
                };
                if end != 0 {
                    x2 = x2.get(..len).unwrap().to_string();
                    x4 = x4.get(..len).unwrap().to_string();
                };
                FQRead(x.0, x2, x.2, x4)
            })
            .collect();
    }

    // Filtering Ns
    chunked_reads.retain(|x| x.1.matches('N').count() <= ns);

    // Filtering Qs
    chunked_reads.retain(|x| {
        let length_cutoff = (x.3.len() as f32 * limit) as usize;
        x.3.as_bytes().iter().filter(|&x| x <= &qv).count() < length_cutoff
    });

    sender.send(chunked_reads).unwrap();
}

pub fn se_writer_worker(recv: Receiver<Vec<FQRead>>, cleanq: String, dedup: bool) {
    let mut clean_file = helper::write_file(cleanq.as_str());
    let mut dups: HashSet<u64> = HashSet::new();

    for r in recv {
        for l in r {
            if dedup {
                let hash = calculate_hash(&l.1);
                if dups.contains(&hash) {
                    continue;
                } else {
                    dups.insert(hash);
                }
            }
            writeln!(clean_file, "{}", l.0).unwrap();
            writeln!(clean_file, "{}", l.1).unwrap();
            writeln!(clean_file, "{}", l.2).unwrap();
            writeln!(clean_file, "{}", l.3).unwrap();
        }
    }
}

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
