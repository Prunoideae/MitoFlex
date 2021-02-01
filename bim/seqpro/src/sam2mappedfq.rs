use super::helper;
use clap::ArgMatches;
use std::collections::HashMap;
use std::path::Path;
use std::string::String;

use std::io::{BufRead, Write};

//Read input as SAM, output {seqname}.fq to directory
pub fn pipeline(arg: ArgMatches) {
    let mut writers: HashMap<String, Box<dyn Write>> = HashMap::new();
    let input = helper::get_input(arg.value_of("input"));
    let output = Path::new(arg.value_of("output").unwrap());

    for l in input
        .lines()
        .map(|x| x.unwrap())
        .filter(|x| !x.starts_with('@'))
    {
        let ls = l.split_whitespace().collect::<Vec<_>>();

        let writer = writers
            .entry(ls.get(2).unwrap().to_string())
            .or_insert_with(|| {
                helper::get_output(Some(
                    output
                        .join(format!("{}.fq", ls.get(2).unwrap()))
                        .to_str()
                        .unwrap(),
                ))
            });

        writeln!(writer, "@{}", ls.get(0).unwrap()).unwrap();
        writeln!(writer, "{}", ls.get(9).unwrap()).unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "{}", ls.get(10).unwrap()).unwrap();
    }
}
