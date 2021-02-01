extern crate clap;
use super::helper;

use clap::ArgMatches;

use std::io::{BufRead, Write};

// Takes input as SAM, output to fq1 and fq2
pub fn pipeline(arg: ArgMatches) {
    let input = helper::get_input(arg.value_of("input"));

    let fq1 = helper::get_output(arg.value_of("fq1"));

    if arg.is_present("fq2") {
        let fq2 = helper::get_output(arg.value_of("fq2"));
        sam2fq2(input, fq1, fq2);
    } else {
        sam2fqi(input, fq1);
    };
}

pub fn sam2fq2(input: Box<dyn BufRead>, mut fq1: Box<dyn Write>, mut fq2: Box<dyn Write>) {
    for l in input
        .lines()
        .map(|x| x.unwrap())
        .filter(|x| !x.starts_with('@'))
    {
        let ls = l.split_whitespace().collect::<Vec<_>>();
        if ls.get(3).unwrap().parse::<usize>().unwrap() & 0x40 != 0 {
            writeln!(fq1, "@{}", ls.get(0).unwrap()).unwrap();
            writeln!(fq1, "{}", ls.get(9).unwrap()).unwrap();
            writeln!(fq1, "+").unwrap();
            writeln!(fq1, "{}", ls.get(10).unwrap()).unwrap();
        } else {
            writeln!(fq2, "@{}", ls.get(0).unwrap()).unwrap();
            writeln!(fq2, "{}", ls.get(9).unwrap()).unwrap();
            writeln!(fq2, "+").unwrap();
            writeln!(fq2, "{}", ls.get(10).unwrap()).unwrap();
        };
    }
}

pub fn sam2fqi(input: Box<dyn BufRead>, mut fq1: Box<dyn Write>) {
    for l in input
        .lines()
        .map(|x| x.unwrap())
        .filter(|x| !x.starts_with('@'))
    {
        let ls = l.split_whitespace().collect::<Vec<_>>();
        writeln!(fq1, "@{}", ls.get(0).unwrap()).unwrap();
        writeln!(fq1, "{}", ls.get(9).unwrap()).unwrap();
        writeln!(fq1, "+").unwrap();
        writeln!(fq1, "{}", ls.get(10).unwrap()).unwrap();
    }
}
