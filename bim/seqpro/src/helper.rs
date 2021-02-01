extern crate flate2;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use std::ffi::OsStr;
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub fn get_input(input: Option<&str>) -> Box<dyn BufRead> {
    if let Some(file_name) = input {
        let path = Path::new(file_name);
        let file = match File::open(&path) {
            Err(_) => panic!("Cannot open file {}!", path.display()),
            Ok(f) => f,
        };

        if path.extension() == Some(OsStr::new("gz")) {
            Box::new(BufReader::with_capacity(128 * 1024, GzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(128 * 1024, file))
        }
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, stdin()))
    }
}

pub fn get_output(output: Option<&str>) -> Box<dyn Write> {
    if let Some(file_name) = output {
        let path = Path::new(file_name);
        let file = match File::create(&path) {
            Err(_) => panic!("Cannot open file {}!", path.display()),
            Ok(f) => f,
        };

        if path.extension() == Some(OsStr::new("gz")) {
            Box::new(BufWriter::with_capacity(
                128 * 1024,
                GzEncoder::new(file, Compression::default()),
            ))
        } else {
            Box::new(BufWriter::with_capacity(128 * 1024, file))
        }
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, stdout()))
    }
}
