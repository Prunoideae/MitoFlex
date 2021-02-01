extern crate clap;

use super::helper;
use clap::ArgMatches;
use itertools::Itertools;
use std::io::BufRead;

// Takes fq1, fq2 as input, and write output to output.
pub fn pipeline(args: ArgMatches) {
    let fq1 = helper::get_input(args.value_of("fq1"));
    let fq2 = helper::get_input(args.value_of("fq2"));
    let mut out = helper::get_output(args.value_of("output"));

    for ((a1, b1), (a2, b2), (a3, b3), (a4, b4)) in fq1.lines().zip(fq2.lines()).tuples() {
        writeln!(out, "{}", a1.unwrap()).unwrap();
        writeln!(out, "{}", a2.unwrap()).unwrap();
        writeln!(out, "{}", a3.unwrap()).unwrap();
        writeln!(out, "{}", a4.unwrap()).unwrap();

        writeln!(out, "{}", b1.unwrap()).unwrap();
        writeln!(out, "{}", b2.unwrap()).unwrap();
        writeln!(out, "{}", b3.unwrap()).unwrap();
        writeln!(out, "{}", b4.unwrap()).unwrap();
    }
}
