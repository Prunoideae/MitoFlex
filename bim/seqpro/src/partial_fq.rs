extern crate clap;
extern crate itertools;

use super::helper;
use clap::ArgMatches;
use itertools::Itertools;
use std::collections::HashSet;
use std::io::{stderr, BufRead, Write};

// Read from input, and extract sequences that are in fq1, but not in input.
pub fn pipeline(arg: ArgMatches) {
    let input = helper::get_input(arg.value_of("input"));
    let fq1 = helper::get_input(arg.value_of("fq1"));

    let mut output = helper::get_output(arg.value_of("output"));
    let mut hashset: HashSet<String> = HashSet::new();

    let mut logger = stderr();

    writeln!(logger, "Loading sequence indexes from input...").unwrap();

    for (index, l1) in input.lines().map(|x| x.unwrap()).enumerate() {
        let ls = l1.split_whitespace().collect::<Vec<_>>();
        if index % 100000 == 0 {
            writeln!(logger, "Loaded {} sequences.", index).unwrap();
        }
        hashset.insert("@".to_string() + ls.get(0).unwrap());
    }

    writeln!(logger, "Load finished, total {} sequences.", hashset.len()).unwrap();

    for (l1, l2, l3, l4) in fq1.lines().map(|x| x.unwrap()).tuples() {
        let ls = l1.split_whitespace().collect::<Vec<_>>();

        if !hashset.contains(&ls.get(0).unwrap().to_string()) {
            writeln!(output, "{}", l1).unwrap();
            writeln!(output, "{}", l2).unwrap();
            writeln!(output, "{}", l3).unwrap();
            writeln!(output, "{}", l4).unwrap();
        }
    }
}
