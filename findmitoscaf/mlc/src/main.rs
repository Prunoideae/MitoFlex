extern crate bcmp;
extern crate bio;
extern crate clap;
extern crate itertools;

use bcmp::{longest_common_substring, AlgoSpec, Match};
use bio::io::fasta::{Reader, Record, Writer};
use clap::{App, Arg};
use itertools::Itertools;
use std::cmp::{max, min};
use std::path::Path;

fn main() {
    let matches = App::new("MLC")
        .author("Junyu Li, <2018301050@szu.edu.cn>")
        .about("Merge some candidate sequences by doing LCS search.")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .help("Input fasta file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .help("Output fasta file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("overlap")
                .short("l")
                .long("overlap")
                .help("Overlapped length to merge")
                .takes_value(true)
                .default_value("10"),
        )
        .get_matches();

    let mut records = Reader::from_file(Path::new(matches.value_of("input").unwrap()))
        .unwrap()
        .records()
        .map(|x| x.unwrap())
        .collect::<Vec<_>>();

    let length = matches
        .value_of("overlap")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    if length < 10 {
        panic!("Min overlap length should be 10 or greater!")
    }

    let mut index = 0;

    loop {
        let mut matched = Match {
            first_pos: 0,
            second_pos: 0,
            length: 0,
        };
        let mut index_i = 0;
        let mut index_j = 0;
        let mut found = false;

        for ((ii, i), (ij, j)) in (&records).iter().enumerate().tuple_combinations() {
            matched = longest_common_substring(i.seq(), j.seq(), AlgoSpec::HashMatch(10));
            if matched.length >= length && (matched.first_pos == 0 || matched.second_pos == 0) {
                println!(
                    "{} {} {} {} {}",
                    ii, ij, matched.first_pos, matched.second_pos, matched.length
                );
                index_i = ii;
                index_j = ij;
                found = true;
                break;
            }
        }

        if !found {
            break;
        }

        if matched.length == records[index_i].seq().len() {
            records.remove(index_i);
        } else if matched.length == records[index_j].seq().len() {
            records.remove(index_j);
        } else {
            let merged_slice = if matched.first_pos > matched.second_pos {
                [
                    &(records[index_i].seq()[0..(matched.first_pos + matched.length)]),
                    &(records[index_j].seq()[(matched.second_pos + matched.length)..]),
                ]
                .concat()
            } else {
                [
                    &(records[index_j].seq()[0..(matched.second_pos + matched.length)]),
                    &(records[index_i].seq()[(matched.first_pos + matched.length)..]),
                ]
                .concat()
            };

            let merged_record = Record::with_attrs(
                &(format!("M{}", index).to_string()),
                Some("flag=1 multi=32767"),
                &merged_slice[..],
            );

            index += 1;

            let (r_i, r_j) = (min(index_i, index_j), max(index_i, index_j) - 1);
            records.remove(r_i);
            records.remove(r_j);

            records.insert(records.len(), merged_record);
        }
    }

    let mut output = Writer::to_file(Path::new(matches.value_of("output").unwrap())).unwrap();
    for r in records {
        output.write_record(&r).ok();
    }

    output.flush().ok();
}
