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

PYBIND11_MODULE(model, m) {

    //PYBIND11_NUMPY_DTYPE(SphCoefs, LMAX, NMAX, NUMT, t, coefs);

    py::class_<MWLMC>(m, "MWLMC")
        .def(py::init<>())
        .def("mwhalo_fields", py::overload_cast<double,double,double,double,int,bool>(&MWLMC::mwhalo_fields),R"pbdoc(
                Return all fields for the MW halo (in the frame of the MW).

                Parameters
                ----------
                t : float
                x, y, z : float
                mwhharmonicflag : int = 127

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

        .def("return_mw_weights", py::overload_cast<double,double,double>(&MWLMC::get_mw_function_weights), R"pbdoc(
             Return function weights for the Milky Way halo.

             Parameters
             ----------
             x,y,z : floats

             Returns
             -------
             potential,fx,fy,fz : array-like float, shaped ((l+1)^2,n)
             )pbdoc",
             py::arg("x"),
             py::arg("y"),
             py::arg("z"))

        .def("return_lmc_weights", py::overload_cast<double,double,double>(&MWLMC::get_lmc_function_weights), R"pbdoc(
             Return function weights for the LMC.

             Parameters
             ----------
             x,y,z : floats

             Returns
             -------
             potential,fx,fy,fz : array-like float, shaped ((l+1)^2,n)
             )pbdoc",
             py::arg("x"),
             py::arg("y"),
             py::arg("z"))

        .def("return_disc_weights", py::overload_cast<double,double,double>(&MWLMC::get_disc_function_weights), R"pbdoc(
             Return function weights for the LMC.

             Parameters
             ----------
             x,y,z: floats

             Returns
             -------
             potential,fx,fy,fz,potential,fx,fy,fz : array-like float, shaped (m+1,n)
             (first set is cosine, second set is sine)
             )pbdoc",
             py::arg("x"),
             py::arg("y"),
             py::arg("z"))

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
             py::arg("discframe")       = false,
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
             py::arg("discframe")       = false,
             py::arg("verbose")         = false)

        .def("print_orbit", &MWLMC::print_orbit, "print an orbit array")

        // no specific resets exposed as of now. could be exposed if this is a common use case.
        //.def("reset_mw_coefficients", &MWLMC::reset_mw_coefficients, "reset MW coefficients")
        .def("reset_all_coefficients", &MWLMC::reset_all_coefficients, "reset all coefficients")



        .def("return_mw_coefficients", &MWLMC::return_mw_coefficients, "return MW halo coefficients")
        .def("return_lmc_coefficients", &MWLMC::return_lmc_coefficients, "return LMC coefficients")
        .def("return_disc_coefficients", &MWLMC::return_disc_coefficients, "return MW disc coefficients")


        // these could possibly be extended to take their own time arrays (which would mean interpolation is possible, finer changes, etc)
        .def("install_mw_coefficients", py::overload_cast<vector<MatrixXd>>(&MWLMC::install_mw_coefficients), "install MW halo coefficients")
        .def("install_lmc_coefficients", py::overload_cast<vector<MatrixXd>>(&MWLMC::install_lmc_coefficients), "install LMC coefficients")
        .def("install_disc_coefficients", py::overload_cast<vector<MatrixXd>,vector<MatrixXd>>(&MWLMC::install_disc_coefficients), "install MW disc coefficients");

}
