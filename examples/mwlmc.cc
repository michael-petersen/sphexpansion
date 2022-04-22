/*

*/

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// local to this build
#include "fullintegrate.h"

namespace py = pybind11;

PYBIND11_MODULE(mwlmc, m) {

    py::class_<MWLMC>(m, "MWLMC")
        .def(py::init<>())
        .def("mwhalo_fields", &MWLMC::mwhalo_fields)
        .def("lmc_fields", &MWLMC::lmc_fields)
        .def("mwd_fields", &MWLMC::mwd_fields)
        .def("all_forces", &MWLMC::all_forces)
        .def("orbit", &MWLMC::orbit)
        .def("print_orbit", &MWLMC::print_orbit);

}
