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
        .def("mwhalo_fields", &MWLMC::mwhalo_fields
             //py::arg("mwhharmonicflag") = 127
            )
        .def("lmc_fields", &MWLMC::lmc_fields)
        .def("mwd_fields", &MWLMC::mwd_fields)
        .def("all_forces", &MWLMC::all_forces)
        .def("get_expansion_centres_physical", &MWLMC::get_expansion_centres_physical)
        .def("get_expansion_centres_virial", &MWLMC::get_expansion_centres_virial)
        .def("orbit", &MWLMC::orbit,
             py::arg("xinit"), py::arg("vinit"),
             py::arg("nint"), py::arg("dt"),
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("fixedtime") = false)
        .def("print_orbit", &MWLMC::print_orbit);

}
