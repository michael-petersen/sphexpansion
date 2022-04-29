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
        .def("mwhalo_fields", &MWLMC::mwhalo_fields,
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("globalframe")     = false,
             py::arg("mwhharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("lmc_fields", &MWLMC::lmc_fields,
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("globalframe")     = false,
             py::arg("lmcharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("mwd_fields", &MWLMC::mwd_fields,
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("globalframe")     = false,
             py::arg("mwdharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("all_forces", &MWLMC::all_forces,
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("globalframe")     = true,
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("get_lmc_trajectory", &MWLMC::get_lmc_trajectory,
             py::arg("dt")     = native_timestep)

        .def("orbit", &MWLMC::orbit,
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("tbegin")          = -2.5,
             py::arg("tend")            = 0.0,
             py::arg("dt")              = 0.002,
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("discframe§")       = true)

        .def("mworbit", &MWLMC::mworbit,
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("tbegin")          = -2.5,
             py::arg("tend")            = 0.0,
             py::arg("dt")              = 0.002)

        .def("rewind", &MWLMC::rewind,
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("dt")              = 0.002,
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("rewindtime")      = 2.5,
             py::arg("discframe")       = true)


        .def("print_orbit", &MWLMC::print_orbit);

}
