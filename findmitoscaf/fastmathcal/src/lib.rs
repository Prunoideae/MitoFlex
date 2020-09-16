#[macro_use]
extern crate pyo3;
use pyo3::prelude::*;

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
    let qs = qs - 1;
    let ss = ss - 1;
    let mut ss = ss;
    let mut se = se;

    if alen >= que || alen >= sub {
        return Ok(true);
    }

    if ss > se {
        let tmp = sub - ss;
        ss = sub - se;
        se = tmp;
    }

    let l = sub + if qs > ss { qe - se } else { se - qe };
    if l > max_length {
        return Ok(false);
    }

    return Ok(l > sub && l > que);
}

#[pymodule]
fn libfastmathcal(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(merge_calculation))?;
    Ok(())
}
