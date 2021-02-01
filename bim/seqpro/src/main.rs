extern crate clap;

mod helper;
mod interleave_fastq;
mod sam2fq;
mod sam2mappedfq;
mod partial_fq;

use clap::{App, Arg};

fn main() {
    let matches = App::new("seqpro")
        .version("v1.0.0")
        .about("Process some sequences in certain ways.")
        .author("Henry Lee")
        .arg(
            Arg::with_name("command")
                .takes_value(true)
                .required(true)
                .display_order(0)
                .possible_values(&["interleave-fastq", "sam2fq", "sam2mappedfq", "partial-fq"])
                .help("What's gonna do with the parameters"),
        )
        .arg(
            Arg::with_name("fq1")
                .short("1")
                .long("fastq1")
                .takes_value(true)
                .help("Fastq file 1")
                .required_ifs(&[
                    ("command", "interleave-fastq"), 
                    ("command", "sam2fq"),
                    ("command","partial_fq"),
                    ]),
        )
        .arg(
            Arg::with_name("fq2")
                .short("2")
                .long("fastq2")
                .takes_value(true)
                .help("Fastq file 2")
                .required_ifs(&[("command", "interleave-fastq")]),
        )
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .takes_value(true)
                .help("Input from, leave out to stdin"),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .takes_value(true)
                .help("Output to (or prefix), leave out to stdout, sometimes handled by other arguments.")
                .required_ifs(&[("command", "sam2mappedfq")]),
                
        )
        .get_matches();

    match matches.value_of("command").unwrap() {
        "interleave-fastq" => interleave_fastq::pipeline(matches),
        "sam2fq" => sam2fq::pipeline(matches),
        "sam2mappedfq" => sam2mappedfq::pipeline(matches),
        "partial-fq"=>partial_fq::pipeline(matches),
        _ => panic!(),
    }
}
