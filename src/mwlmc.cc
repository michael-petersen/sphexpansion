/*


*/

// generic header includes: requires both pybind11 and eigen
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// eigen includes
#include <Eigen/StdVector>
#include <Eigen/Dense>
using Eigen::MatrixXd;

// local to this build
#include "fullintegrate.h"

namespace py = pybind11;

PYBIND11_MODULE(mwlmc, m) {
    // py::options options;
    // options.disable_function_signatures();

    py::class_<MWLMC>(m, "MWLMC")
        .def(py::init<>())
        .def("mwhalo_fields", py::overload_cast<double,double,double,double,int,bool>(&MWLMC::mwhalo_fields),R"pbdoc(
                Return all fields for the MW halo (in the frame of the MW).

                Parameters
                ----------
                t : float
                x, y, z : float
                lmcharmonicflag : int = 127

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
             py::arg("mwhharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("mwhalo_fields", py::overload_cast<double,std::vector<double>,std::vector<double>,std::vector<double>,int,bool>(&MWLMC::mwhalo_fields),R"pbdoc(
                Return all fields for the MW halo (in the frame of the MW halo).

                Parameters
                ----------
                t : float
                x, y, z : array-like float
                lmcharmonicflag : int = 127

                Returns
                -------
                fx, fy, fz : array-like float
                density : array-like float
                potential : array-like float
            )pbdoc",
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("mwhharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("lmc_fields", py::overload_cast<double,std::vector<double>,std::vector<double>,std::vector<double>,int,bool>(&MWLMC::lmc_fields),R"pbdoc(
                Return all fields for the LMC halo (in the frame of the LMC).

                Parameters
                ----------
                t : float
                x, y, z : array-like float
                lmcharmonicflag : int = 127

                Returns
                -------
                fx, fy, fz : array-like float
                density : array-like float
                potential : array-like float
            )pbdoc",
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("lmcharmonicflag") = 127,
             py::arg("verbose")         = false)

       .def("lmc_fields", py::overload_cast<double,double,double,double,int,bool>(&MWLMC::lmc_fields),R"pbdoc(
                Return all fields for the LMC halo (in the frame of the LMC).

                Parameters
                ----------
                t : float
                x, y, z : float
                lmcharmonicflag : int = 127

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
            //py::arg("globalframe")     = false,
            py::arg("lmcharmonicflag") = 127,
            py::arg("verbose")         = false)

        .def("mwd_fields", py::overload_cast<double,std::vector<double>,std::vector<double>,std::vector<double>,int,bool>(&MWLMC::mwd_fields),R"pbdoc(
                Return all fields for the MW disc (in the frame of the MW disc).

                Parameters
                ----------
                t : float
                x, y, z : array-like float
                lmcharmonicflag : int = 127

                Returns
                -------
                fx, fy, fz : array-like float
                density : array-like float
                potential : array-like float
            )pbdoc",
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             //py::arg("globalframe")     = false,
             py::arg("mwdharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("mwd_fields", py::overload_cast<double,double,double,double,int,bool>(&MWLMC::mwd_fields),R"pbdoc(
                Return all fields for the MW disc (in the frame of the MW disc).

                Parameters
                ----------
                t : float
                x, y, z : float
                lmcharmonicflag : int = 127

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
             py::arg("mwdharmonicflag") = 127,
             py::arg("verbose")         = false)

        .def("all_forces", py::overload_cast<double,double,double,double,int,int,int,bool>(&MWLMC::all_forces), R"pbdoc(
                Return total forces (in the frame of the MW disc).

                Parameters
                ----------
                t : float
                x, y, z : float
                mwhharmonicflag : int = 127
                mwdharmonicflag : int = 127
                lmcharmonicflag : int = 127

                Returns
                -------
                fx, fy, fz : float
            )pbdoc",
             py::arg("t"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"),
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("verbose")         = false)

       .def("all_forces", py::overload_cast<double,std::vector<double>,std::vector<double>,std::vector<double>,int,int,int,bool>(&MWLMC::all_forces),R"pbdoc(
                Return total forces (in the frame of the MW disc).

                Parameters
                ----------
                t : float
                x, y, z : array-like float
                mwhharmonicflag : int = 127
                mwdharmonicflag : int = 127
                lmcharmonicflag : int = 127

                Returns
                -------
                fx, fy, fz : array-like float
            )pbdoc",
            py::arg("t"),
            py::arg("x"),
            py::arg("y"),
            py::arg("z"),
            //py::arg("globalframe")     = true,
            py::arg("mwhharmonicflag") = 127,
            py::arg("mwdharmonicflag") = 127,
            py::arg("lmcharmonicflag") = 127,
            py::arg("verbose")         = false)

        .def("expansion_centres", py::overload_cast<double,bool>(&MWLMC::get_expansion_centres_physical),R"pbdoc(
                 Get the expansion centres at a specific moment in simulation time

                 Parameters
                 ----------
                 t : float

                 Returns
                 -------
                 tx,ty,tz,mwhx,mwhy,mwhz,lmcx,lmcy,lmcz,mwdx,mwdy,mwdz : array of expansion centres
             )pbdoc",
             py::arg("t"),
             py::arg("verbose")         = false)

       .def("expansion_centre_velocities", py::overload_cast<double,bool>(&MWLMC::get_expansion_centre_velocities_physical),R"pbdoc(
                Get the expansion centres at a specific moment in simulation time

                Parameters
                ----------
                t : float

                Returns
                -------
                tx,ty,tz,mwhx,mwhy,mwhz,lmcx,lmcy,lmcz,mwdx,mwdy,mwdz : array of expansion centres
            )pbdoc",
            py::arg("t"),
            py::arg("verbose")         = false)

        .def("get_lmc_trajectory", &MWLMC::get_lmc_trajectory,R"pbdoc(
             Get the LMC trajectory (relative to the MW disc centre).
                Parameters
                ----------
                rewindtime : float = 2.5
                dt : float = native_timestep

                Returns
                -------
                t : array-like float
                x, y, z : array-like float
                )pbdoc",
             py::arg("rewindtime") = 2.5,
             py::arg("dt")     = native_timestep)

        .def("mworbit",  py::overload_cast<vector<double>,vector<double>,double,double,double>(&MWLMC::mworbit), R"pbdoc(
             Compute an orbit integration using only the initial Milky Way potential.

             Parameters
             ----------
             xinit : array-like float
             vinit : array-like float
             tbegin : float = -2.5
             tend : float = 0.
             dt : float = 0.002

             Returns
             -------
             x, y, z : array-like float
             vx, vy, vz :array-like float
             fx, fy, fz :array-like float
             t :array-like float
             )pbdoc",
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("tbegin")          = -2.5,
             py::arg("tend")            = 0.0,
             py::arg("dt")              = 0.002)

        .def("mworbit", py::overload_cast<MatrixXd,MatrixXd,double,double,double>(&MWLMC::mworbit), R"pbdoc(
             Compute an orbit integration using only the initial Milky Way potential.

             Parameters
             ----------
             xinit : array-like float
             vinit : array-like float
             tbegin : float = -2.5
             tend : float = 0.
             dt : float = 0.002

             Returns
             -------
             x, y, z : array-like float
             vx, vy, vz : array-like float
             fx, fy, fz : array-like float
             t : array-like float
             )pbdoc",
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("tbegin")          = -2.5,
             py::arg("tend")            = 0.0,
             py::arg("dt")              = 0.002)

        .def("rewind", py::overload_cast<vector<double>,vector<double>,double,int,int,int,double,bool,bool>(&MWLMC::rewind),R"pbdoc(
             Compute an orbit rewind.

             Parameters
             ----------
             xinit : array-like float
             vinit : array-like float
             dt : float = 0.002
             mwhharmonicflag : int = 127
             mwdharmonicflag : int = 127
             lmcharmonicflag : int = 127
             rewindtime : float = 2.5

             Returns
             -------
             x, y, z : array-like float
             vx, vy, vz : array-like float
             fx, fy, fz : array-like float
             t : array-like float
             )pbdoc",
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("dt")              = 0.002,
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("rewindtime")      = 2.5,
             py::arg("discframe")       = true,
             py::arg("verbose")         = false)

        .def("rewind", py::overload_cast<MatrixXd,MatrixXd,double,int,int,int,double,bool,bool>(&MWLMC::rewind), R"pbdoc(
             Compute an orbit rewind.

             Parameters
             ----------
             xinit : array-like float
             vinit : array-like float
             dt : float = 0.002
             mwhharmonicflag : int = 127
             mwdharmonicflag : int = 127
             lmcharmonicflag : int = 127
             rewindtime : float = 2.5

             Returns
             -------
             x, y, z : array-like float
             vx, vy, vz : array-like float
             fx, fy, fz : array-like float
             t : array-like float
             )pbdoc",
             py::arg("xinit"),
             py::arg("vinit"),
             py::arg("dt")              = 0.002,
             py::arg("mwhharmonicflag") = 127,
             py::arg("mwdharmonicflag") = 127,
             py::arg("lmcharmonicflag") = 127,
             py::arg("rewindtime")      = 2.5,
             py::arg("discframe")       = true,
             py::arg("verbose")         = false)

        .def("print_orbit", &MWLMC::print_orbit, "print an orbit array");

}
