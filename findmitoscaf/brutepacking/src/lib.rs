#[macro_use]
extern crate cpython;

use cpython::{PyResult, Python, PyList};

py_module_initializer!(libbrutepacking, initlibbrutepacking, PyInit_libbrutepacking, |py, m| {
    m.add(
        py,
        "__doc__",
        "A rust solution for brute-forcing set packing problems.",
    )?;
    m.add(py, "solution", py_fn!(py, solution_py(sets:PyList)))?;
    Ok(())
});

fn solution_py(_: Python, sets:PyList) -> PyResult<PyList> {
    Ok(sets)
}
