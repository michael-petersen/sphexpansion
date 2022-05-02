/*


*/

// generic header includes: requires both pybind11 and eigen
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// local to this build
#include "fullintegrate.h"

namespace py = pybind11;

PYBIND11_MODULE(mwlmc, m) {
    // py::options options;
    // options.disable_function_signatures();

    py::class_<MWLMC>(m, "MWLMC")
        .def(py::init<>())
        .def("mwhalo_fields", &MWLMC::mwhalo_fields,
             "Return all fields for the MW halo (in the frame of the MW halo).",
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("globalframe")     = false,
             py::arg("mwhharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("lmc_fields", &MWLMC::lmc_fields, R"pbdoc(
                Return all fields for the LMC halo (in the frame of the LMC).

                Parameters
                ----------
                t : float
                x, y, z : float
                globalframe : bool = False
                lmcharmonicflag : int = 127
                verbose : bool = False

                Returns
                
                -------
                fx, fy, fz : float
                density : float
                potential : float

            )pbdoc",
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("globalframe")     = false,
             py::arg("lmcharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("mwd_fields", &MWLMC::mwd_fields,
             "Return all fields for the MW disc.",
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("globalframe")     = false,
             py::arg("mwdharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("all_forces", &MWLMC::all_forces,
             "Return total forces (in the frame of the MW disc).",
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
             "Get the LMC trajectory (relative to the MW disc centre).",
             py::arg("dt")     = native_timestep)

        .def("mworbit", &MWLMC::mworbit,
             "compute an orbit integration using only the initial \
              Milky Way potential.",
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("tbegin")          = -2.5,
             py::arg("tend")            = 0.0,
             py::arg("dt")              = 0.002)

        .def("mworbit_parallel", &MWLMC::mworbit_parallel,
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("tbegin")          = -2.5,
             py::arg("tend")            = 0.0,
             py::arg("dt")              = 0.002)

        .def("rewind", &MWLMC::rewind, R"pbdoc(
                Compute an orbit rewind in all three components.
                
                Parameters
                ----------
                xinit : array-like
                vinit : array-like
                dt : 
                
                Returns
                -------
            )pbdoc",
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("dt")              = 0.002,
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("rewindtime")      = 2.5,
             py::arg("discframe")       = true)

        .def("print_orbit", &MWLMC::print_orbit, "print an orbit array");

}
