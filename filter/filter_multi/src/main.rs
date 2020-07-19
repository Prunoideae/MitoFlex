extern crate itertools;

mod app;
mod helper;
mod pe;
mod se;

use std::sync::mpsc::sync_channel;
use std::thread::spawn;
use std::thread::JoinHandle;

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

    let dedup: bool = matches.is_present("deduplication");
    let trunc: bool = matches.is_present("truncate");

    if pe {
        filter_pe(
            fastq1, fastq2, cleanq1, cleanq2, start, end, ns, quality, limit, dedup, trim, trunc,
        )
    } else {
        filter_se(
            fastq1, cleanq1, start, end, ns, quality, limit, trim, trunc, dedup,
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
) {
    let fastq1 = String::from(fastq1);
    let cleanq1 = String::from(cleanq1);

    let mut handles = vec![];

    // Reading file, output : recv_w
    let (send_w, recv_w) = sync_channel::<se::FQRead>(1000);
    handles.push(spawn(move || se::se_read_file(fastq1, send_w, trim)));

    // Trimming file, output : recv_t (if trimmed) / recv_w (if not trimmed)
    let recv_t = if start != 0 || end != 0 {
        let (send_t, recv_t) = sync_channel::<se::FQRead>(1000);
        handles.push(spawn(move || se::se_trim_read(recv_w, start, end, send_t)));
        recv_t
    } else {
        recv_w
    };

    let recv_d = if !trunc {
        // Filtering ns, output : recv_n
        let (send_n, recv_n) = sync_channel::<se::FQRead>(1000);
        handles.push(spawn(move || se::se_filter_n(recv_t, ns, send_n)));

        // Filtering qualities, output : recv_q
        let (send_q, recv_q) = sync_channel::<se::FQRead>(1000);
        handles.push(spawn(move || {
            se::se_filter_q(recv_n, quality, limit, send_q)
        }));

        // Filtering duplicates, output : recv_d (if deduped) / recv_q (if not)
        if dedup {
            let (send_d, recv_d) = sync_channel::<se::FQRead>(1000);
            handles.push(spawn(move || se::se_filter_d(recv_q, send_d)));
            recv_d
        } else {
            recv_q
        }
    } else {
        recv_t
    };

    // Writing to file
    handles.push(spawn(move || se::se_write_file(recv_d, cleanq1)));

    //Join handles
    for i in handles {
        i.join().unwrap();
    }
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
) {
    let fastq1 = String::from(fastq1);
    let fastq2 = String::from(fastq2);

    let cleanq1 = String::from(cleanq1);
    let cleanq2 = String::from(cleanq2);

    let mut handles: Vec<JoinHandle<()>> = vec![];

    // Reading files, output : recv_w
    let (send_w, recv_w) = sync_channel::<pe::FQRead2>(1000);
    handles.push(spawn(move || {
        pe::pe_read_file(fastq1, fastq2, send_w, trim)
    }));

    // Trimming files, output : recv_t (if trimmed) / recv_w (if not)
    let recv_t = if start != 0 || end != 0 {
        let (send_t, recv_t) = sync_channel::<pe::FQRead2>(1000);
        handles.push(spawn(move || pe::pe_trim_read(recv_w, start, end, send_t)));
        recv_t
    } else {
        recv_w
    };

    let recv_d = if !trunc {
        //Filtering Ns, output : recv_n
        let (send_n, recv_n) = sync_channel::<pe::FQRead2>(1000);
        handles.push(spawn(move || pe::pe_filter_n(recv_t, ns, send_n)));

        //Filtering Qs, output : recv_q
        let (send_q, recv_q) = sync_channel::<pe::FQRead2>(1000);
        handles.push(spawn(move || {
            pe::pe_filter_q(recv_n, quality, limit, send_q)
        }));

        //Filtering duplicates, output : recv_d
        if dedup {
            let (send_d, recv_d) = sync_channel::<pe::FQRead2>(1000);
            handles.push(spawn(move || pe::pe_filter_d(recv_q, send_d)));
            recv_d
        } else {
            recv_q
        }
    } else {
        recv_t
    };

    //Writing to files
    handles.push(spawn(move || pe::pe_write_file(recv_d, cleanq1, cleanq2)));

    for i in handles {
        i.join().unwrap();
    }
}
