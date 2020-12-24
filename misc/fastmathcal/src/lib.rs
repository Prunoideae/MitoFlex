#[macro_use]
extern crate pyo3;

use bio::alphabets::dna;
use bio::io::fasta;
use ndarray::prelude::*;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[inline(always)]
fn merge_calculation_internal(
    que: usize,
    sub: usize,
    alen: usize,
    qs: usize,
    qe: usize,
    ss: usize,
    se: usize,
    max_length: usize,
) -> bool {
    let qs = qs - 1;
    let ss = ss - 1;
    let mut ss = ss;
    let mut se = se;

    if alen >= que || alen >= sub {
        return true;
    }

    if ss > se {
        let tmp = sub - ss;
        ss = sub - se;
        se = tmp;
    }

    let l = sub + if qs > ss { qe - se } else { se - qe };
    if l > max_length {
        return false;
    }

    return l > sub && l > que;
}

#[pyfunction]
fn merge_calculation(
    que: usize,
    sub: usize,
    alen: usize,
    qs: usize,
    qe: usize,
    ss: usize,
    se: usize,
    max_length: usize,
) -> PyResult<bool> {
    Ok(merge_calculation_internal(
        que, sub, alen, qs, qe, ss, se, max_length,
    ))
}

#[pyfunction]
fn wash_merge_blast(
    blast_file: &str,
    contigs_file: &str,
    search_range: isize,
    overlapped_len: usize,
    max_length: usize,
) -> PyResult<()> {
    let result_file = format!("{}.filtered", blast_file);
    let blast_path = Path::new(blast_file);
    let blast_in = BufReader::with_capacity(1024 * 128, File::open(blast_path)?);
    let out_path = Path::new(&result_file);
    let mut blast_out = BufWriter::with_capacity(1024 * 128, File::create(out_path)?);

    let mut seqslen: HashMap<String, usize> = HashMap::new();
    for seq in fasta::Reader::from_file(Path::new(contigs_file))
        .unwrap()
        .records()
        .map(|x| x.unwrap())
    {
        seqslen.insert(seq.id().to_string(), seq.seq().len());
    }

    let mut seqseq: HashSet<String> = HashSet::new();

    for l in blast_in.lines().map(|x| x.unwrap()) {
        let row = l.split('\t').collect::<Vec<_>>();

        let query = row.get(0).unwrap();
        let subject = row.get(1).unwrap();
        let alen = row.get(3).unwrap().parse::<usize>()?;
        let qs = row.get(6).unwrap().parse::<usize>()?;
        let qe = row.get(7).unwrap().parse::<usize>()?;
        let ss = row.get(8).unwrap().parse::<usize>()?;
        let se = row.get(9).unwrap().parse::<usize>()?;

        if query == subject {
            continue;
        }

        if alen < overlapped_len {
            continue;
        }

        if seqseq.contains(&format!("{}{}", subject, query)) {
            continue;
        }
        seqseq.insert(format!("{}{}", query, subject));

        let que = *seqslen.get(&query.to_string()).unwrap_or(&0);
        let sub = *seqslen.get(&subject.to_string()).unwrap_or(&0);
        if que == 0 || sub == 0 {
            continue;
        }

        if alen < que && alen < sub && search_range >= 0 {
            if (ss > search_range as usize && sub - se > search_range as usize)
                || (qs > search_range as usize && que - qe > search_range as usize)
            {
                continue;
            }
        }

        if !merge_calculation_internal(que, sub, alen, qs, qe, ss, se, max_length) {
            continue;
        }

        writeln!(blast_out, "{}", l)?;
    }
    Ok(())
}

#[pyfunction]
fn merge_overlaps(
    blast_file: &str,
    contigs_file: &str,
    out_file: &str,
    sindex: usize,
) -> PyResult<usize> {
    let out_path = Path::new(out_file);
    let mut contig_out = fasta::Writer::new(BufWriter::with_capacity(
        128 * 1024,
        File::create(out_path)?,
    ));
    let mut sindex = sindex;

    let mut seqs: HashMap<String, fasta::Record> = HashMap::new();
    fasta::Reader::from_file(Path::new(contigs_file))
        .unwrap()
        .records()
        .map(|x| x.unwrap())
        .for_each(|x| {
            seqs.insert(x.id().to_string(), x);
        });

    let mut blast_records = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(Path::new(blast_file))
        .unwrap()
        .records()
        .map(|x| x.unwrap())
        .collect::<Vec<_>>();

    while !blast_records.is_empty() {
        let overlapped = blast_records.pop().unwrap();
        let que = overlapped.get(0).unwrap();
        let subj = overlapped.get(1).unwrap();
        let sque = seqs.remove(&que.to_string()).unwrap();
        let mut ssub = seqs.remove(&subj.to_string()).unwrap();
        let alen = overlapped.get(3).unwrap().parse::<usize>().unwrap();
        let qs: u32 = overlapped.get(6).unwrap().parse::<u32>().unwrap() - 1;
        let qe: u32 = overlapped.get(7).unwrap().parse().unwrap();

        let mut ss: u32 = overlapped.get(8).unwrap().parse::<u32>().unwrap() - 1;
        let mut se: u32 = overlapped.get(9).unwrap().parse().unwrap();

        if ss > se {
            let tmp: u32 = ssub.seq().len() as u32 - ss;
            ss = ssub.seq().len() as u32 - (se - 1);
            se = tmp;
            ssub = fasta::Record::with_attrs(ssub.id(), ssub.desc(), &dna::revcomp(ssub.seq())[..]);
            ssub.seq().clone_from(&&dna::revcomp(ssub.seq())[..]);
        }

        let new_seq = if alen >= sque.seq().len() {
            ssub
        } else if alen >= ssub.seq().len() {
            sque
        } else {
            let sec_conc = &if qs > ss {
                let mut new_seq = sque.seq()[..qe as usize]
                    .iter()
                    .map(|x| *x)
                    .collect::<Vec<_>>();
                new_seq.extend_from_slice(&ssub.seq()[se as usize..]);
                new_seq
            } else {
                let mut new_seq = ssub.seq()[..se as usize]
                    .iter()
                    .map(|x| *x)
                    .collect::<Vec<_>>();
                new_seq.extend_from_slice(&sque.seq()[qe as usize..]);
                new_seq
            }[..];

            fasta::Record::with_attrs(
                &format!("M{}", sindex),
                Some(&format!("flag=1 multi=32767 len={}", sec_conc.len())),
                sec_conc,
            )
        };

        contig_out.write_record(&new_seq)?;
        blast_records.retain(|x| {
            x.get(0).unwrap() != que
                && x.get(1).unwrap() != que
                && x.get(0).unwrap() != subj
                && x.get(1).unwrap() != subj
        });
        sindex += 1;
    }

    seqs.values()
        .for_each(|x| contig_out.write_record(x).unwrap());
    Ok(sindex)
}

#[pyfunction]
fn seq_overlap(seq1: &str, seq2: &str) -> PyResult<(usize, usize, usize)> {
    let n = seq1.len();
    let m = seq2.len();
    let seq1 = seq1.as_bytes();
    let seq2 = seq2.as_bytes();

    let mut dp = Array2::<usize>::zeros((n, m));
    let mut maxi = 0;
    let mut maxj = 0;
    let mut maxv = 0;

    for i in 0..n {
        for j in 0..m {
            let curr_val = if seq1.get(i) == seq2.get(j) {
                dp.get((i - 1, j - 1)).unwrap_or(&0) + 1
            } else {
                0
            };
            *dp.get_mut((i, j)).unwrap() = curr_val;
            if curr_val > maxv {
                maxi = i;
                maxj = j;
                maxv = curr_val;
            }
        }
    }

    Ok((maxi + 1 - maxv, maxj + 1 - maxv, maxv))
}

#[pymodule]
fn libfastmathcal(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(wash_merge_blast, m)?)?;
    m.add_function(wrap_pyfunction!(merge_overlaps, m)?)?;
    m.add_function(wrap_pyfunction!(merge_calculation, m)?)?;
    m.add_function(wrap_pyfunction!(seq_overlap, m)?)?;
    Ok(())
}
