extern crate clap;
use clap::{App, Arg, ArgMatches};
use std::cmp::PartialOrd;
use std::fmt::Display;
use std::fs::metadata;
use std::ops::Bound;
use std::ops::RangeBounds;
use std::path::Path;
use std::str::FromStr;

fn format_bounds<T, D>(r: D) -> String
where
    T: PartialOrd + FromStr + Display,
    D: RangeBounds<T>,
{
    format!(
        "{}..{}",
        match r.start_bound() {
            Bound::Unbounded => "".to_string(),
            Bound::Included(n) => n.to_string(),
            Bound::Excluded(n) => n.to_string(),
        },
        match r.end_bound() {
            Bound::Unbounded => "".to_string(),
            Bound::Included(n) => n.to_string(),
            Bound::Excluded(n) => n.to_string(),
        }
    )
}

fn ranged_validator<T, D>(x: String, r: D) -> Result<(), String>
where
    T: PartialOrd + FromStr + Display,
    D: RangeBounds<T>,
{
    if let Ok(n) = x.parse::<T>() {
        if r.contains(&n) {
            return Ok(());
        }
    }
    Err(format!("{} is not in {}", x, format_bounds(r)))
}

pub fn get_app() -> ArgMatches<'static> {
    App::new("Fastq Filter")
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
                .takes_value(true)
                .validator(|x| {
                    if let Ok(f) = metadata(&x) {
                        if f.is_file() {
                            return Ok(());
                        }
                    }
                    Err(format!("{} is not a valid path for reading file!", x))
                }),
        )
        .arg(
            Arg::with_name("fastq2")
                .short("2")
                .long("fastq2")
                .value_name("FASTQ2")
                .help("Input raw data fastq file 2")
                .takes_value(true)
                .validator(|x| {
                    if let Ok(f) = metadata(&x) {
                        if f.is_file() {
                            return Ok(());
                        }
                    }
                    Err(format!("{} is not a valid path for reading file!", x))
                }),
        )
        .arg(
            Arg::with_name("cleanq1")
                .short("3")
                .long("cleanq1")
                .value_name("CLEANQ1")
                .help("Output clean fastq file 1")
                .required(true)
                .takes_value(true)
                .validator(|x| {
                    let dirname = match Path::new(&x).parent() {
                        Some(n) => n,
                        None => return Err(String::from("Not a valid path")),
                    };

                    if let Ok(f) = metadata(&dirname) {
                        if f.is_dir() {
                            return Ok(());
                        }
                    }
                    Err(format!("{} is not a valid path for reading file!", x))
                }),
        )
        .arg(
            Arg::with_name("cleanq2")
                .short("4")
                .long("cleanq2")
                .value_name("CLEANQ2")
                .help("Output clean fastq file 2")
                .requires("fastq2")
                .takes_value(true)
                .validator(|x| {
                    let dirname = match Path::new(&x).parent() {
                        Some(n) => n,
                        None => return Err(String::from("Not a valid path")),
                    };

                    if let Ok(f) = metadata(&dirname) {
                        if f.is_dir() {
                            return Ok(());
                        }
                    }
                    Err(format!("{} is not a valid path for reading file!", x))
                }),
        )
        .arg(
            Arg::with_name("start")
                .short("s")
                .long("start")
                .value_name("INT")
                .help("Cut sequence's start")
                .takes_value(true)
                .default_value("0")
                .validator(|x| ranged_validator(x, 0..)),
        )
        .arg(
            Arg::with_name("end")
                .short("e")
                .long("end")
                .value_name("INT")
                .help("Cut suquences's end")
                .takes_value(true)
                .default_value("0")
                .validator(|x| ranged_validator(x, 0..)),
        )
        .arg(
            Arg::with_name("quality")
                .short("q")
                .long("quality")
                .value_name("INT")
                .help("Quality under this will be considered as bad base")
                .takes_value(true)
                .default_value("55")
                .validator(|x| ranged_validator(x, 0..)),
        )
        .arg(
            Arg::with_name("limit")
                .short("l")
                .long("limit")
                .value_name("FLOAT")
                .help("Sequences will be filtered out if bad bases > limit * length")
                .takes_value(true)
                .default_value("0.2")
                .validator(|x| ranged_validator(x, 0.0..1.0)),
        )
        .arg(
            Arg::with_name("nvalues")
                .short("n")
                .long("nvalues")
                .value_name("INT")
                .help("Sequences having Ns more than this will be filtered out")
                .takes_value(true)
                .default_value("10")
                .validator(|x| ranged_validator(x, 0..)),
        )
        .arg(
            Arg::with_name("trimming")
                .short("t")
                .long("trim")
                .value_name("INT")
                .help("Only this bases of sequences will be filtered out")
                .takes_value(true)
                .default_value("0")
                .validator(|x| ranged_validator(x, 0..)),
        )
        .arg(
            Arg::with_name("deduplication")
                .short("d")
                .long("deduplication")
                .help("Filter out duplicated sequences")
                .requires("fastq2"),
        )
        .arg(
            Arg::with_name("truncate")
                .long("truncate_only")
                .help("Only truncates the file, no filtering."),
        )
        .get_matches()
}
