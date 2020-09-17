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
            Arg::with_name("depth")
                .short("d")
                //.required(true)
                .takes_value(true)
                .help("min depth")
                .value_name("INT"),
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
        .arg(
            Arg::with_name("deny_number")
                .short("m")
                .takes_value(true)
                .help("Use it to take only x sequence with highest depth.")
                .value_name("INT")
                .conflicts_with("depth"),
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
    let depth: i32 = app.value_of("depth").unwrap().parse().unwrap();

    let mut count = 0;

    if !app.is_present("deny_number") {
        for (l1, l2) in infile.lines().tuples() {
            let title = l1.unwrap();
            let seq = l2.unwrap();

            if !title.starts_with(">") {
                continue;
            }

            if depth != 0 {
                let seq_depth = title.split_whitespace().collect::<Vec<&str>>()[2]
                    .to_string()
                    .split('=')
                    .collect::<Vec<&str>>()[1]
                    .parse::<f32>()
                    .unwrap();

                if depth as f32 > seq_depth {
                    continue;
                }
            }

            let length = seq.len() - 1;
            if length < min || length > max {
                continue;
            }
            writeln!(outfile, "{}", title).unwrap();
            writeln!(outfile, "{}", seq).unwrap();
            count += 1;
        }
    } else {
        let max_count = app
            .value_of("deny_number")
            .unwrap()
            .parse::<usize>()
            .unwrap();

        let seqs = infile
            .lines()
            .map(|x| x.unwrap())
            .tuples::<(_, _)>()
            .filter(|(_, b)| b.len() <= max && b.len() >= min)
            .collect::<Vec<_>>();
        seqs.sort_by_cached_key(|(a, _)| {
            a.split_whitespace().collect::<Vec<_>>()[2]
                .split('=')
                .collect::<Vec<_>>()[1]
                .parse::<f32>()
                .unwrap();
        });
        seqs.take(max_count).for_each(|(t, s)| {
            count += 1;
            writeln!(outfile, "{}", t).unwrap();
            writeln!(outfile, "{}", s).unwrap();
        });
    }

    println!("{}", count);
}
