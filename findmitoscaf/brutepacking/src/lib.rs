/*
    lib brutepacking
    =========
    Copyright (c) 2019-2020 Li Junyu <2018301050@szu.edu.cn>.
    This file is part of MitoFlex.
    MitoFlex is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    MitoFlex is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with MitoFlex.  If not, see <http://www.gnu.org/licenses/>.
*/

#[macro_use]
extern crate cpython;
extern crate itertools;

use cpython::{PyResult, Python};
use std::clone::Clone;
use std::cmp::{max, min};

#[derive(Clone)]
struct Unit {
    sets: Vec<(u32, u32)>,
    keys: Vec<Vec<(u32, u32)>>,
}

impl Unit {
    fn new(init: Vec<(u32, u32)>, key: Vec<Vec<(u32, u32)>>) -> Unit {
        Unit {
            sets: init,
            keys: key,
        }
    }

    fn push(&mut self, key: Vec<(u32, u32)>) {
        if self.sets.is_empty() {
            self.sets = key.clone();
        } else {
            let mut new_sets: Vec<(u32, u32)> = Vec::new();
            for (i, r) in key.iter().zip(&self.sets) {
                if i.0 == 0 && i.1 == 0 {
                    new_sets.push((r.0, r.1))
                } else {
                    let start = min(i.0, r.0);
                    let end = max(i.1, r.1);
                    new_sets.push((start, end));
                }
            }
            self.sets = new_sets;
        }
        self.keys.push(key);
    }

    fn check(&self, grate: &Vec<(u32, u32)>) -> bool {
        for (i, r) in grate.iter().zip(&self.sets) {
            if i.0 == 0 && i.1 == 0 {
                continue;
            }
            if overlapped((i.0, i.1), (r.0, r.1)) {
                return false;
            }
        }
        true
    }

    fn value(&self) -> u32 {
        let mut total = 0;
        for r in &self.sets {
            total += r.1 - r.0;
        }
        total
    }
}

py_module_initializer!(brute, initbrute, PyInit_brute, |py, m| {
    m.add(
        py,
        "__doc__",
        "A rust solution for brute-forcing set packing problems.",
    )?;
    m.add(py, "__version__", "1.3")?;
    m.add(
        py,
        "solution",
        py_fn!(py, solution_py(sets: Vec<Vec<(u32, u32)>>)),
    )?;
    Ok(())
});

fn overlapped(r1: (u32, u32), r2: (u32, u32)) -> bool {
    max(r1.1, r2.1) - min(r1.0, r2.0) < ((r1.1 - r1.0) + (r2.1 - r2.0))
}

/*
    Use a brute-force method to break through everything. Proved time complexity
    for this is n to 2^n-1, where the worst situtaion is all the sets are disjoint.
    Since low number of PCGs required and sequence output by genewise is always slow.
*/
fn solution(sets: Vec<Vec<(u32, u32)>>) -> Vec<Vec<(u32, u32)>> {
    let mut pool: Vec<Unit> = Vec::new();
    //Push a empty seed into it
    pool.push(Unit::new(Vec::new(), Vec::new()));

    for set in &sets {
        let mut next_gen: Vec<Unit> = Vec::new();
        for seed in &pool {
            if seed.check(set) {
                let mut son = seed.clone();
                son.push(set.clone());
                next_gen.push(son);
            }
        }
        pool = [pool, next_gen].concat();
    }

    let mut result: Unit = Unit::new(Vec::new(), Vec::new());
    for r in pool {
        println!("{:?}", r.keys);
        if r.value() > result.value() {
            result = r;
        }
    }
    result.keys
}

#[test]
fn test_solution() {
    let a = vec![
        vec![(1, 3), (5, 6)],
        vec![(2, 3), (0, 0)],
        vec![(4, 7), (9, 10)],
    ];
    let b = vec![vec![(2, 3), (0, 0)], vec![(4, 7), (9, 10)]];
    assert_eq!(solution(a), b);
}

fn solution_py(_: Python, sets: Vec<Vec<(u32, u32)>>) -> PyResult<Vec<Vec<(u32, u32)>>> {
    Ok(solution(sets))
}
