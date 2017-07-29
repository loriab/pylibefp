/*
    pylibefp/core.cc: Main binding of libefp with pybind11

    Copyright (c) 2017 The Psi4 Developers

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/


#include "efp.h"
#include "private.h"
//#include "pybind11.h"

#include <pybind11/pybind11.h>

//int add(int i, int j) {
//    return i + j;
//}
//
//int subtract(int i, int j) {
//    return i - j;
//}

std::string wrapped_efp_banner(efp* efp) {
    std::string str;
    str = std::string(efp_banner());
    return str;
}

void wrapped_efp_print_banner(efp* efp) {
    printf("%s", efp_banner());
}

int wrapped_get_frag_atom_count(efp* efp, int frag_idx) {
    enum efp_result res;
    size_t n=0;

    res = efp_get_frag_atom_count(efp, frag_idx, &n);
    //if ((res = efp_get_frag_atom_count(efp_, frag_idx, &n)))
    //    throw PsiException("EFP::get_frag_atom_count(): " +
    //        std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return n;
}




namespace py = pybind11;

PYBIND11_PLUGIN(pylibefp) {
    py::module m("pylibefp", "Python wrapping of Parallel implementation of the Effective Fragment Potential (EFP) method");

//    m.def("add", &add, R"pbdoc(
//        Add two numbers
//        Some other explanation about the add function.
//    )pbdoc");
//
//    m.def("subtract", &subtract, R"pbdoc(
//        Subtract two numbers
//        Some other explanation about the subtract function.
//    )pbdoc");

    m.attr("__version__") = py::str("1.1");

    py::enum_<efp_result>(m, "efp_result", "Result of a libefp operation")
        .value("EFP_RESULT_SUCCESS", EFP_RESULT_SUCCESS)                      /* Operation was successful. */
        .value("EFP_RESULT_FATAL", EFP_RESULT_FATAL)                          /* Fatal error has occurred. */
        .value("EFP_RESULT_NO_MEMORY", EFP_RESULT_NO_MEMORY)                  /* Insufficient memory. */
        .value("EFP_RESULT_FILE_NOT_FOUND", EFP_RESULT_FILE_NOT_FOUND)        /* File not found. */
        .value("EFP_RESULT_SYNTAX_ERROR", EFP_RESULT_SYNTAX_ERROR)            /* Syntax error. */
        .value("EFP_RESULT_UNKNOWN_FRAGMENT", EFP_RESULT_UNKNOWN_FRAGMENT)    /* Unknown EFP fragment. */
        .value("EFP_RESULT_POL_NOT_CONVERGED", EFP_RESULT_POL_NOT_CONVERGED)  /* Polarization SCF procedure did not converge. */
        .export_values();

    py::enum_<efp_term>(m, "efp_term", "Flags to specify EFP energy terms")
        .value("EFP_TERM_ELEC", EFP_TERM_ELEC)        /* EFP/EFP electrostatics. */
        .value("EFP_TERM_POL", EFP_TERM_POL)          /* EFP/EFP polarization. */
        .value("EFP_TERM_DISP", EFP_TERM_DISP)        /* EFP/EFP dispersion. */
        .value("EFP_TERM_XR", EFP_TERM_XR)            /* EFP/EFP exchange repulsion. */
        .value("EFP_TERM_CHTR", EFP_TERM_CHTR)        /* EFP/EFP charge transfer, reserved for future use. */
        .value("EFP_TERM_AI_ELEC", EFP_TERM_AI_ELEC)  /* Ab initio/EFP electrostatics. */
        .value("EFP_TERM_AI_POL", EFP_TERM_AI_POL)    /* Ab initio/EFP polarization. */
        .value("EFP_TERM_AI_DISP", EFP_TERM_AI_DISP)  /* Ab initio/EFP dispersion, reserved for future use. */
        .value("EFP_TERM_AI_XR", EFP_TERM_AI_XR)      /* Ab initio/EFP exchange repulsion, reserved for future use. */
        .value("EFP_TERM_AI_CHTR", EFP_TERM_AI_CHTR)  /* Ab initio/EFP charge transfer, reserved for future use. */
        .export_values();

    py::enum_<efp_disp_damp>(m, "efp_disp_damp", "Fragment-fragment dispersion damping type")
        .value("EFP_DISP_DAMP_OVERLAP", EFP_DISP_DAMP_OVERLAP)  /* Overlap-based damping (default). */
        .value("EFP_DISP_DAMP_TT", EFP_DISP_DAMP_TT)            /* Tang-Toennies damping. */
        .value("EFP_DISP_DAMP_OFF", EFP_DISP_DAMP_OFF)          /* No dispersion damping. */
        .export_values();

    py::enum_<efp_elec_damp>(m, "efp_elec_damp", "Fragment-fragment electrostatic damping type")
        .value("EFP_ELEC_DAMP_SCREEN", EFP_ELEC_DAMP_SCREEN)    /* SCREEN-controlled damping (default). */
        .value("EFP_ELEC_DAMP_OVERLAP", EFP_ELEC_DAMP_OVERLAP)  /* Overlap-based damping. */
        .value("EFP_ELEC_DAMP_OFF", EFP_ELEC_DAMP_OFF)          /* No electrostatic damping. */
        .export_values();

    py::enum_<efp_pol_damp>(m, "efp_pol_damp", "Fragment-fragment polarization damping type")
        .value("EFP_POL_DAMP_TT", EFP_POL_DAMP_TT)    /* Tang-Toennies like damping (default). */
        .value("EFP_POL_DAMP_OFF", EFP_POL_DAMP_OFF)  /* No polarization damping. */
        .export_values();

    py::enum_<efp_coord_type>(m, "efp_coord_type", "Describes the way fragment coordinates are specified")
        .value("EFP_COORD_TYPE_XYZABC", EFP_COORD_TYPE_XYZABC)  /* Coordinates of center of mass of a fragment and Euler angles. */
        .value("EFP_COORD_TYPE_POINTS", EFP_COORD_TYPE_POINTS)  /* Coordinates of three points belonging to a fragment. */
        .value("EFP_COORD_TYPE_ROTMAT", EFP_COORD_TYPE_ROTMAT)  /* Coordinates of fragment center of mass and its rotation matrix. */
        .export_values();


    py::enum_<efp_pol_driver>(m, "efp_pol_driver", "Driver used for solving polarization equations")
        .value("EFP_POL_DRIVER_ITERATIVE", EFP_POL_DRIVER_ITERATIVE)  /* Iterative solution of polarization equations. */
        .value("EFP_POL_DRIVER_DIRECT", EFP_POL_DRIVER_DIRECT)        /* Direct solution of polarization equations. */
        .export_values();

    py::class_<efp_opts>(m, "efp_opts", "Options controlling EFP computation")
        .def(py::init());
        /* Specifies which energy terms to compute. */
        /* Dispersion damping type (see #efp_disp_damp). */
        /* Electrostatic damping type (see #efp_elec_damp). */
        /* Polarization damping type (see #efp_pol_damp). */
        /* Driver used to find polarization induced dipoles. */
        /* Enable periodic boundary conditions if nonzero. */
        /* Enable fragment-fragment interaction cutoff if nonzero. */
        /* Cutoff distance for fragment-fragment interactions. */
        // unsigned terms;                      /* Specifies which energy terms to compute. */
        // enum efp_disp_damp disp_damp;        /* Dispersion damping type (see #efp_disp_damp). */
        // enum efp_elec_damp elec_damp;        /* Electrostatic damping type (see #efp_elec_damp). */
        // enum efp_pol_damp pol_damp;          /* Polarization damping type (see #efp_pol_damp). */
        // enum efp_pol_driver pol_driver;      /* Driver used to find polarization induced dipoles. */
        // int enable_pbc;                      /* Enable periodic boundary conditions if nonzero. */
        // int enable_cutoff;                   /* Enable fragment-fragment interaction cutoff if nonzero. */
        // double swf_cutoff;                   /* Cutoff distance for fragment-fragment interactions. */

    py::class_<efp_energy>(m, "efp_energy", "EFP energy terms")
        .def(py::init());
        /* EFP/EFP electrostatic energy. */
        /* Charge penetration energy from overlap-based electrostatic damping. Zero if overlap-based damping is turned off. */
        /* Interaction energy of EFP electrostatics with point charges. */
        /* All polarization energy goes here. Polarization is computed self-consistently so it can't be separated into EFP/EFP and AI/EFP parts. */
        /* EFP/EFP dispersion energy. */
        /* AI/EFP dispersion energy. */
        /* EFP/EFP exchange-repulsion energy. */
        /* Sum of all the above energy terms. */
        // double electrostatic;
        // double charge_penetration;
        // double electrostatic_point_charges;
        // double polarization;
        // double dispersion;
        // double ai_dispersion;
        // double exchange_repulsion;
        // double total;

    py::class_<efp_atom>(m, "efp_atom", "EFP atom info")
        .def(py::init());
        // char label[32];   /**< Atom label. */
        // double x;         /**< X coordinate of atom position. */
        // double y;         /**< Y coordinate of atom position. */
        // double z;         /**< Z coordinate of atom position. */
        // double mass;      /**< Atom mass. */
        // double znuc;      /**< Nuclear charge. */

    py::class_<efp>(m, "efp", "Main libefp opaque structure")
        .def(py::init())
        //.def("result", &efp_result, "Callback function which is called by libefp to obtain electric field *arg2* in the number *arg0* specified points at positions *arg1* with initialized user data *arg3*")
        .def("banner", &efp_banner, "Gets a human readable banner string with information about the library")
        .def("print_banner", &efp_print_banner, "Prints libefp banner to stdout")
        .def("efp_create", &efp_create, "Creates a new efp object")
        .def("opts_default", &efp_opts_default, "Gets default values of simulation options returns in *arg0*")
        .def("set_error_log", &efp_set_error_log, "Sets the error log callback function")
        .def("set_opts", &efp_set_opts, "Set computation options to *arg0*")
        .def("get_opts", &efp_get_opts, "Gets currently set computation options returns in *arg0*")
        .def("add_potential", &efp_add_potential, "Adds EFP potential from a file *arg0*")
        .def("add_fragment", &efp_add_fragment, "Adds a new fragment *arg0* to the EFP subsystem")
        .def("prepare", &efp_prepare, "Prepares the calculation")
        .def("skip_fragments", &efp_skip_fragments, "Skip interactions between the fragments *arg0* and *arg1* inclusive if *arg2*")
        .def("set_electron_density_field_fn", &efp_set_electron_density_field_fn, "Sets the callback function which computes electric field from electrons in ab initio subsystem to *arg0*")
        .def("set_electron_density_field_user_data", &efp_set_electron_density_field_user_data, "Sets user data *arg0* to be passed to ::efp_electron_density_field_fn")
        .def("set_point_charges", &efp_set_point_charges, "Setup *arg0* arbitrary point charges of magnitude *arg1* at locations *arg2* interacting with EFP subsystem")
        .def("get_point_charge_count", &efp_get_point_charge_count, "Gets the number of currently set point charges and return it in *arg0*")
        .def("get_point_charge_values", &efp_get_point_charge_values, "Gets values of currently set point charges and returns them in *arg1*")
        .def("set_point_charge_values", &efp_set_point_charge_values, "Sets values of point charges *arg0*")
        .def("get_point_charge_coordinates", &efp_get_point_charge_coordinates, "Gets coordinates of currently set point charges and returns them in *arg1*")
        .def("set_point_charge_coordinates", &efp_set_point_charge_coordinates, "Sets coordinates *arg0* of point charges")
        .def("get_point_charge_gradient", &efp_get_point_charge_gradient, "Gets gradient on point charges from EFP subsystem and returns them in *arg1*")
        .def("set_coordinates", &efp_set_coordinates, "Update positions and orientations of all fragments with types in array *arg0* and returns them in *arg1*")
        .def("set_frag_coordinates", &efp_set_frag_coordinates, "Updates position and orientation of the specified effective 0-indexed fragment *arg0* of type *arg1* and returns it in *arg2*")
        .def("get_coordinates", &efp_get_coordinates, "Gets center of mass positions and Euler angles of the effective fragments and returns it in *arg0*")
        .def("get_frag_xyzabc", &efp_get_frag_xyzabc, "Gets center of mass position and Euler angles on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("set_periodic_box", &efp_set_periodic_box, "Sets up periodic box size of *arg0* by *arg1* by *arg2*")
        .def("get_stress_tensor", &efp_get_stress_tensor, "Gets the stress tensor and returns it in *arg0*")
        .def("get_ai_screen", &efp_get_ai_screen, "Gets the ab initio screening parameters on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("set_orbital_energies", &efp_set_orbital_energies, "Sets ab initio orbital energies to *efp0* number core orbitals, *efp1* number active orbitals, *efp2* number virtual orbitals, *efp3* array of orbital energies")
        .def("set_dipole_integrals", &efp_set_dipole_integrals, "Sets ab initio dipole integrals to  *efp0* number core orbitals, *efp1* number active orbitals, *efp2* number virtual orbitals, *efp3* dipole integral matrices")
        .def("get_wavefunction_dependent_energy", &efp_get_wavefunction_dependent_energy, "Updates wavefunction-dependent energy terms returning in *arg0*")
        .def("compute", &efp_compute, "Perform the EFP computation, doing gradient if *arg0*")
        .def("get_frag_charge", &efp_get_frag_charge, "Gets total charge on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_frag_charge", &efp_get_frag_charge, "Gets total charge on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_frag_multiplicity", &efp_get_frag_multiplicity, "Gets spin multiplicity on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_frag_multiplole_count", &efp_get_frag_multipole_count, "Gets number of electrostatic multipole points on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_multipole_count", &efp_get_multipole_count, "Gets total number of multipoles from EFP electrostatics and returns it in *arg0*")
        .def("get_multipole_coordinates", &efp_get_multipole_coordinates, "Gets coordinates of electrostatics multipoles and returns it in *arg0*")
        .def("get_multipole_values", &efp_get_multipole_values, "Gets electrostatics multipoles from EFP fragments and returns it in *arg0*")
        .def("get_induced_dipole_count", &efp_get_induced_dipole_count, "Gets the number of polarization induced dipoles and returns it in *arg0*")
        .def("get_induced_dipole_coordinates", &efp_get_induced_dipole_coordinates, "Gets coordinates of induced dipoles and returns it in *arg0*")
        .def("get_induced_dipole_values", &efp_get_induced_dipole_values, "Gets values of polarization induced dipoles and returns it in *arg0*")
        .def("get_induced_dipole_conj_values", &efp_get_induced_dipole_conj_values, "Gets values of polarization conjugated induced dipoles and returns it in *arg0*")
        .def("get_lmo_count", &efp_get_lmo_count, "Gets the number of LMOs in a fragment and returns it in *arg0*")
        .def("get_lmo_coordinates", &efp_get_lmo_coordinates, "Gets coordinates of LMO centroids on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_xrfit", &efp_get_xrfit, "Gets parameters of fitted exchange-repulsion on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_energy", &efp_get_energy, "Gets computed energy components and returns it in *arg0*")
        .def("get_gradient", &efp_get_gradient, "Gets computed EFP energy gradient and returns it in *arg0*")
        .def("get_atomic_gradient", &efp_get_atomic_gradient, "Gets computed EFP energy gradient on individual atoms and returns it in *arg0*")
        .def("get_frag_count", &efp_get_frag_count, "Gets the number of fragments in this computation and returns it in *arg0*")
        .def("get_frag_name", &efp_get_frag_name, "Gets the name of the specified 0-indexed effective fragment *arg0* and returns it in *arg2* of length *arg1*")
        .def("get_frag_mass", &efp_get_frag_mass, "Gets total mass on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_frag_inertia", &efp_get_frag_inertia, "Gets fragment principal moments of inertia on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_frag_atom_count", &efp_get_frag_atom_count, "Gets the number of atoms on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_frag_atoms", &efp_get_frag_atoms, "Gets atoms comprising the specified 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_electric_field", &efp_get_electric_field, "Gets electric field for a point on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("torque_to_derivative", &efp_torque_to_derivative, "Convert rigid body torque *arg1* to derivatives *arg2* of energy by Euler angles *arg0*")
        .def("shutdown", &efp_shutdown, "Release all resources used by this EFP")
        .def("result_to_string", &efp_result_to_string, "Result value to be converted to string");

    return m.ptr();
}

        //.def("banner", wrapped_efp_banner, "Gets a human readable banner string with information about the library")
        //.def("print_banner", wrapped_efp_print_banner, "Prints libefp banner to stdout")

        //.def("get_frag_atom_count", wrapped_get_frag_atom_count, "Gets the number of atoms on 0-indexed fragment *arg0* and returns it in *arg1*")

