/*

*/

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "fullintegrate.h"

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE(mwlmc, m) {
//     m.doc() = "pybind11 example plugin"; // optional module docstring

//     m.def("add", &add, "A function that adds two numbers");

    py::class_<MWLMC>(m, "MWLMC")
        .def(py::init<>())
        .def("mwhalo_fields", &MWLMC::mwhalo_fields)
        .def("lmc_fields", &MWLMC::lmc_fields)
        .def("mwd_fields", &MWLMC::mwd_fields)
        .def("all_forces", &MWLMC::all_forces)
        .def("orbit", &MWLMC::orbit)
        .def("print_orbit", &MWLMC::print_orbit);

}
