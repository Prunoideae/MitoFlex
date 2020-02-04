extern crate clap;
extern crate flate2;
extern crate itertools;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub fn read_file(file_name: &str) -> Box<dyn BufRead> {
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
}

pub fn write_file(file_name: &str) -> Box<dyn Write> {
    let path = Path::new(file_name);
    let file = match File::create(&path) {
        Err(_) => panic!("Cannot open file {}", path.display()),
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
}
