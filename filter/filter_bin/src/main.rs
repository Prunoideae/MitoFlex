extern crate clap;
extern crate flate2;
extern crate itertools;

mod helper;

use clap::{App, Arg};
use itertools::Itertools;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;

fn main() {
    let matches = App::new("Fastq Filter")
        .version("0.1")
        .author("Junyu Li, <2018301050@szu.edu.cn>")
        .about("Filter out unqualified fastq sequences")
        .arg(
            Arg::with_name("fastq1")
                .short("1")
                .long("fastq1")
                .value_name("FASTQ1")
                .help("Input raw data fastq file 1")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("fastq2")
                .short("2")
                .long("fastq2")
                .value_name("FASTQ2")
                .help("Input raw data fastq file 2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("cleanq1")
                .short("3")
                .long("cleanq1")
                .value_name("CLEANQ1")
                .help("Output clean fastq file 1")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("cleanq2")
                .short("4")
                .long("cleanq2")
                .value_name("CLEANQ2")
                .help("Output clean fastq file 2")
                .requires("fastq2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("start")
                .short("s")
                .long("start")
                .value_name("INT")
                .help("Cut sequence's start")
                .takes_value(true)
                .default_value("0"),
        )
        .arg(
            Arg::with_name("end")
                .short("e")
                .long("end")
                .value_name("INT")
                .help("Cut suquences's end")
                .takes_value(true)
                .default_value("0"),
        )
        .arg(
            Arg::with_name("quality")
                .short("q")
                .long("quality")
                .value_name("INT")
                .help("Quality under this will be considered as bad base")
                .takes_value(true)
                .default_value("55"),
        )
        .arg(
            Arg::with_name("limit")
                .short("l")
                .long("limit")
                .value_name("FLOAT")
                .help("Sequences will be filtered out if bad bases > limit * length")
                .takes_value(true)
                .default_value("0.2"),
        )
        .arg(
            Arg::with_name("nvalues")
                .short("n")
                .long("nvalues")
                .value_name("INT")
                .help("Sequences having Ns more than this will be filtered out")
                .takes_value(true)
                .default_value("10"),
        )
        .arg(
            Arg::with_name("trimming")
                .short("t")
                .long("trim")
                .value_name("INT")
                .help("Only this bases of sequences will be filtered out")
                .takes_value(true)
                .default_value("0"),
        )
        .arg(
            Arg::with_name("deduplication")
                .short("d")
                .long("deduplication")
                .help("Filter out duplicated sequences")
                .requires("fastq2"),
        )
        .get_matches();

    let fastq1 = matches.value_of("fastq1").unwrap();
    let fastq2 = matches.value_of("fastq2").unwrap_or("");
    let cleanq1 = matches.value_of("cleanq1").unwrap();
    let cleanq2 = matches.value_of("cleanq2").unwrap_or("");

    let start: usize = matches
        .value_of("start")
        .unwrap_or("-1")
        .parse()
        .ok()
        .expect("Cannot parse start position!");
    let end: usize = match matches.value_of("end").unwrap_or("-1").parse() {
        Ok(n) => {
            if start > n {
                panic!("Start position comes after the end!")
            } else {
                n
            }
        }
        Err(_) => panic!("Cannot parse end position!"),
    };
    let quality: u8 = match matches.value_of("quality").unwrap_or("55").parse() {
        Ok(n) => {
            if n <= 0 || n > 100 {
                panic!("Wrong quality number!")
            }
            n
        }
        Err(_) => panic!("Canoot parse quality value!"),
    };
    let limit: f32 = match matches.value_of("limit").unwrap_or("0.2").parse() {
        Ok(n) => {
            if n <= 0.0 || n >= 1.0 {
                panic!("Wrong percentage value!")
            }
            n
        }
        Err(_) => panic!("Cannot parse limit value!"),
    };

    let ns: usize = matches
        .value_of("nvalues")
        .unwrap_or("10")
        .parse()
        .ok()
        .expect("Cannot parse N value!");
    let trim: usize = matches
        .value_of("trimming")
        .unwrap_or("0")
        .parse()
        .ok()
        .expect("Cannot parse a positive int to trimming!");

    let dedup: bool = matches.is_present("deduplication");

    if fastq2.is_empty() {
        filter_se(fastq1, cleanq1, start, end, ns, quality, limit, trim);
    } else {
        filter_pe(
            fastq1, fastq2, cleanq1, cleanq2, start, end, ns, quality, limit, dedup, trim,
        );
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
) {
    let fq1 = helper::read_file(fastq1);
    let fq2 = helper::read_file(fastq2);

    let mut cl1 = helper::write_file(cleanq1);
    let mut cl2 = helper::write_file(cleanq2);

    let mut dup: HashSet<u64> = HashSet::new();

    let len = end - start;

    let mut counts = 0;

    for ((a1, b1), (a2, b2), _, (a4, b4)) in fq1.lines().zip(fq2.lines()).tuples() {
        let head1: String = a1.ok().unwrap().trim().to_string();
        let head2: String = b1.ok().unwrap().trim().to_string();
        let mut seq1: String = a2.ok().unwrap().trim().to_string();
        let mut seq2: String = b2.ok().unwrap().trim().to_string();
        let mut qua1: String = a4.ok().unwrap().trim().to_string();
        let mut qua2: String = b4.ok().unwrap().trim().to_string();

        if start != 0 {
            seq1 = seq1.get(start..).unwrap().to_string();
            seq2 = seq2.get(start..).unwrap().to_string();
            qua1 = qua1.get(start..).unwrap().to_string();
            qua2 = qua2.get(start..).unwrap().to_string();
        }

        if end != 0 {
            seq1 = seq1.get(..len).unwrap().to_string();
            seq2 = seq2.get(..len).unwrap().to_string();
            qua1 = qua1.get(..len).unwrap().to_string();
            qua2 = qua2.get(..len).unwrap().to_string();
        }

        if seq1.matches("N").count() > ns || seq2.matches("N").count() > ns {
            continue;
        }

        let cutoff = (seq1.len() as f32 * limit) as usize;
        if qua1.as_bytes().iter().filter(|&n| n <= &quality).count() >= cutoff
            || qua2.as_bytes().iter().filter(|&n| n <= &quality).count() >= cutoff
        {
            continue;
        }

        if dedup {
            let hash = calculate_hash(&seq1);

            if dup.contains(&hash) {
                continue;
            }
            dup.insert(hash);
        }

        if trim != 0 {
            counts += seq1.len();
            if counts > trim {
                break;
            }
        }
        writeln!(cl1, "{}", head1);
        writeln!(cl1, "{}", seq1);
        writeln!(cl1, "+");
        writeln!(cl1, "{}", qua1);
        writeln!(cl2, "{}", head2);
        writeln!(cl2, "{}", seq2);
        writeln!(cl2, "+");
        writeln!(cl2, "{}", qua2);
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
) {
    let fastq_file = helper::read_file(fastq1);
    let mut clean_file = helper::write_file(cleanq1);
    let len = end - start;
    let mut times = 0;
    for (l1, l2, _, l4) in fastq_file.lines().tuples() {
        let head: String = l1.ok().unwrap().trim().to_string();
        let mut bps: String = l2.ok().unwrap().trim().to_string();
        let mut quas: String = l4.ok().unwrap().trim().to_string();

        if start != 0 {
            bps = bps.get(start..).unwrap().to_string();
            quas = quas.get(start..).unwrap().to_string();
        }
        if end != 0 {
            bps = bps.get(..len).unwrap().to_string();
            quas = quas.get(..len).unwrap().to_string();
        }
        if bps.matches("N").count() > ns {
            continue;
        }

        let cutoff = quas.len() as f32 * limit;
        if quas.as_bytes().iter().filter(|&n| n <= &quality).count() >= cutoff as usize {
            continue;
        }

        if trim != 0 {
            times += bps.len();
            if times > trim {
                break;
            }
        }

        writeln!(clean_file, "{}", head);
        writeln!(clean_file, "{}", bps);
        writeln!(clean_file, "+");
        writeln!(clean_file, "{}", quas);
    }
}

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
