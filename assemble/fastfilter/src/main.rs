extern crate clap;

mod helper;

use clap::{App, Arg};
use itertools::Itertools;
use std::io::prelude::*;

fn main() {
    let app = App::new("Length filter")
        .version("0.1")
        .author("Junyu Li")
        .about("A simple filter for trimming intermediate contigs.")
        .arg(
            Arg::with_name("length")
                .short("l")
                .required(true)
                .require_delimiter(true)
                .takes_value(true)
                .help("required min and max length")
                .value_name("INT,INT"),
        )
        .arg(
            Arg::with_name("input")
                .short("i")
                .required(true)
                .takes_value(true)
                .help("Input file path, accepts only one-line format")
                .value_name("PATH"),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .required(true)
                .takes_value(true)
                .help("Output file path")
                .value_name("PATH"),
        )
        .get_matches();

    let lengths = app
        .values_of("length")
        .unwrap()
        .map(|x| x.parse::<usize>().unwrap())
        .collect::<Vec<_>>();

    if lengths.len() != 2 {
        println!("Input length string not valid, please input INT,INT.")
    }

    let (min, max) = (lengths[0], lengths[1]);

    let infile = helper::read_file(app.value_of("input").unwrap());
    let mut outfile = helper::write_file(app.value_of("output").unwrap());

    for (l1, l2) in infile.lines().tuples() {
        let title = l1.ok().unwrap();
        let seq = l2.ok().unwrap();

        let length = seq.len() - 1;
        if length < min || length > max {
            continue;
        }
        write!(outfile, "{}", title).unwrap();
        write!(outfile, "{}", seq).unwrap();
    }

    println!("Hello, world!");
}
