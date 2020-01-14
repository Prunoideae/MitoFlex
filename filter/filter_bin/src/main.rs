extern crate clap;
extern crate flate2;
extern crate itertools;

mod helper;

use clap::{App, Arg};
use flate2::read::GzDecoder;
use itertools::Itertools;
use std::fs::File;
use std::io;
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
            Arg::with_name("times")
                .short("t")
                .long("times")
                .value_name("INT")
                .help("Only this number of sequences will be filtered out")
                .takes_value(true)
                .default_value("-1"),
        )
        .arg(
            Arg::with_name("deduplication")
                .short("d")
                .long("deduplication")
                .help("Filter out duplicated sequences"),
        )
        .get_matches();

    let fastq1 = matches.value_of("fastq1").unwrap();
    let fastq2 = matches.value_of("fastq2").unwrap_or("");
    let cleanq1 = matches.value_of("cleanq1").unwrap();
    let cleanq2 = matches.value_of("cleanq2").unwrap_or("");

    let start: usize = match matches.value_of("start").unwrap_or("-1").parse() {
        Ok(n) => {
            if n < 0 {
                0
            } else {
                n
            }
        }
        Err(_) => panic!("Cannot parse start position!"),
    };
    let end: usize = match matches.value_of("end").unwrap_or("-1").parse() {
        Ok(n) => {
            if n <= 0 {
                0
            } else {
                if start > n {
                    panic!("Start position comes after the end!")
                } else {
                    n
                }
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
    let times: i32 = matches
        .value_of("times")
        .unwrap_or("-1")
        .parse()
        .ok()
        .expect("Cannot parse argument times");

    let dedup: bool = matches.is_present("deduplication");

    if fastq2.is_empty() {
        filter_se(fastq1, cleanq1, start, end, ns, quality, limit);
    } else {
        filter_pe();
    }
}

fn filter_pe() {}

fn filter_se(
    fastq1: &str,
    cleanq1: &str,
    start: usize,
    end: usize,
    ns: usize,
    quality: u8,
    limit: f32,
) {
    let fastq_file = helper::read_file(fastq1);
    let mut clean_file = helper::write_file(cleanq1);
    let len = end - start;
    for mut lines in &fastq_file.lines().chunks(4) {
        let head: String = lines.next().unwrap().unwrap().trim().to_string();
        let mut bps: String = lines.next().unwrap().unwrap().trim().to_string();
        let plus: String = lines.next().unwrap().unwrap().trim().to_string();
        let mut quas: String = lines.next().unwrap().unwrap().trim().to_string();
        let cutoff = quas.len() as f32 * limit;

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

        if quas.as_bytes().iter().filter(|&n| n <= &quality).count() >= cutoff as usize {
            continue;
        }

        writeln!(clean_file, "{}", head);
        writeln!(clean_file, "{}", bps);
        writeln!(clean_file, "{}", plus);
        writeln!(clean_file, "{}", quas);
    }
}

#[test]
fn s() {
    let a = 2;
    assert_eq!(Some("12"), String::from("1234").get(0..0));
}
