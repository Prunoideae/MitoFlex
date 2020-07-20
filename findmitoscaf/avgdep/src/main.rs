extern crate clap;

use clap::{App, Arg};
use std::collections::HashMap;
use std::fs::File;
use std::io::stdin;
use std::io::stdout;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn get_input(input: Option<&str>) -> Box<dyn BufRead> {
    match input {
        Some(n) => {
            let path = Path::new(n);
            Box::new(BufReader::with_capacity(
                1024 * 128,
                match File::open(&path) {
                    Err(_) => panic!("Cannot open file {}!", path.display()),
                    Ok(n) => n,
                },
            ))
        }
        None => Box::new(BufReader::with_capacity(1024 * 128, stdin())),
    }
}

fn get_output(output: Option<&str>) -> Box<dyn Write> {
    match output {
        Some(n) => {
            let path = Path::new(n);
            Box::new(BufWriter::with_capacity(
                1024 * 128,
                match File::open(&path) {
                    Err(_) => panic!("Cannot open file {}!", path.display()),
                    Ok(n) => n,
                },
            ))
        }
        None => Box::new(BufWriter::with_capacity(1024 * 128, stdout())),
    }
}

struct DepthInfo {
    length: usize,
    sum_depth: usize,
}

impl DepthInfo {
    fn new() -> Self {
        DepthInfo {
            length: 0,
            sum_depth: 0,
        }
    }
}

fn main() {
    let matches = App::new("Avg depth calculator")
        .version("0.1")
        .author("Junyu Li, <2018301050@szu.edu.cn>")
        .about("Calculate depth from samtools depth piping.")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .help("Input data from, leave out to stdin."),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .takes_value(true)
                .help("Output data to, leave out to stdout."),
        )
        .get_matches();

    let fin = get_input(matches.value_of("input"));
    let mut fout = get_output(matches.value_of("output"));

    let mut mapping: HashMap<String, DepthInfo> = HashMap::new();

    for i in fin.lines().map(|x| x.unwrap()) {
        let tup: Vec<_> = i.split_whitespace().collect();

        let mut record = mapping
            .entry(tup.get(0).unwrap().to_string())
            .or_insert(DepthInfo::new());

        record.length += 1;
        record.sum_depth += tup.get(2).unwrap().parse::<usize>().unwrap();
    }
    for (k, v) in mapping {
        writeln!(fout, "{} {}", k, v.sum_depth / v.length).unwrap();
    }
}
