extern crate futures;
extern crate futures_cpupool;
extern crate itertools;
extern crate rayon;

mod app;
mod helper;
mod pe;
mod se;

use rayon::join;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::mpsc::sync_channel;
use std::sync::Arc;

fn main() {
    let matches = app::get_app();

    let fastq1 = matches.value_of("fastq1").unwrap();
    let cleanq1 = matches.value_of("cleanq1").unwrap();

    let pe = matches.is_present("fastq2");
    let fastq2 = matches.value_of("fastq2").unwrap_or_default();
    let cleanq2 = matches.value_of("cleanq2").unwrap_or_default();

    let start: usize = matches.value_of("start").unwrap().parse().unwrap();
    let end: usize = matches.value_of("end").unwrap().parse().unwrap();
    let quality: u8 = matches.value_of("quality").unwrap().parse().unwrap();
    let limit: f32 = matches.value_of("limit").unwrap().parse().unwrap();
    let ns: usize = matches.value_of("nvalues").unwrap().parse().unwrap();
    let trim: usize = matches.value_of("trimming").unwrap().parse().unwrap();
    let threads: usize = matches.value_of("threads").unwrap().parse().unwrap();

    let dedup: bool = matches.is_present("deduplication");
    let trunc: bool = matches.is_present("truncate");

    if pe {
        filter_pe(
            fastq1, fastq2, cleanq1, cleanq2, start, end, ns, quality, limit, dedup, trim, trunc,
            threads,
        )
    } else {
        filter_se(
            fastq1, cleanq1, start, end, ns, quality, limit, trim, trunc, dedup, threads,
        )
    }
}

fn filter_se(
    fastq1: &str,
    cleanq1: &str,
    start: usize,
    end: usize,
    ns: usize,
    quality: u8,
    limit: f32,
    trim: usize,
    trunc: bool,
    dedup: bool,
    threads: usize,
) {
    if trunc {
        todo!();
    }
    let fastq1 = String::from(fastq1);
    let cleanq1 = String::from(cleanq1);
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let (send, recv) = sync_channel::<Vec<se::FQRead>>(10000);
    let send = Arc::new(send);
    join(
        move || se::se_reader_worker(fastq1, trim, pool, send, ns, quality, limit, start, end),
        move || se::se_writer_worker(recv, cleanq1, dedup),
    );
}

fn filter_pe(
    fastq1: &str,
    fastq2: &str,
    cleanq1: &str,
    cleanq2: &str,
    start: usize,
    end: usize,
    ns: usize,
    quality: u8,
    limit: f32,
    dedup: bool,
    trim: usize,
    trunc: bool,
    threads: usize,
) {
}
