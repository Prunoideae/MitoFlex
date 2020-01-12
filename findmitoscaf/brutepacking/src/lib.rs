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

use cpython::{PyResult, Python};
use std::clone::Clone;

#[derive(Clone)]
struct Unit {
    sets: Vec<u32>,
    key: u32,
}

impl Unit {
    fn new(init: Vec<u32>, key: u32) -> Unit {
        Unit {
            sets: init,
            key: key,
        }
    }

    fn push(&mut self, grate: u32) {
        self.sets.push(grate);
        self.key |= grate;
    }

    fn check(&self, grate: u32) -> bool {
        self.key & grate == 0
    }
}

py_module_initializer!(brute, initbrute, PyInit_brute, |py, m| {
    m.add(
        py,
        "__doc__",
        "A rust solution for brute-forcing set packing problems.",
    )?;
    m.add(py, "solution", py_fn!(py, solution_py(sets: Vec<u32>)))?;
    Ok(())
});

/*
    Use a brute-force method to break through everything. Proved time complexity
    for this is n to 2^n-1, where the worst situtaion is all the sets are disjoint.
    Since low number of PCGs required and sequence output by genewise is always slow.
*/
fn solution(sets: Vec<u32>) -> Vec<u32> {
    let mut pool: Vec<Unit> = Vec::new();
    pool.push(Unit::new(Vec::new(), 0));
    for i in sets {
        let mut pool_append: Vec<Unit> = Vec::new();
        for u in &pool {
            if u.check(i) {
                let mut new_key: Unit = u.clone();
                new_key.push(i);
                pool_append.push(new_key);
            }
        }
        pool = [pool, pool_append].concat();
    }

    let mut result: Unit = Unit::new(vec![], 0);
    for r in pool {
        if r.key.count_ones() > result.key.count_ones() {
            result = r;
        }
    }
    result.sets
}

#[test]
fn test_solution() {
    assert_eq!(solution(vec![3, 24, 7]), vec![24, 7]);
}

fn solution_py(_: Python, sets: Vec<u32>) -> PyResult<Vec<u32>> {
    Ok(solution(sets))
}
