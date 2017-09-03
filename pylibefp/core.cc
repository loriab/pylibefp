/*
    pylibefp/core.cc: Main binding of libefp with pybind11

    Copyright (c) 2017 The Psi4 Developers

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/


#include <efp.h>
#include <libefp/private.h>

#include <pybind11/pybind11.h>

#include <sstream>
#include <algorithm>

namespace py = pybind11;

class libefpException : public std::exception {
public:
    explicit libefpException(const char * m) : message{m} {}
    virtual const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
};




std::string _efp_banner(efp* efp) {
    std::string str;
    str = std::string(efp_banner());
    return str;
}

py::tuple _efp_get_frag_name(efp* efp, size_t frag_idx) {
    enum efp_result res;
    char buffer[80];

    res = efp_get_frag_name(efp, frag_idx, 80, buffer);
    std::string fname = std::string(buffer);

    py::tuple rets = py::make_tuple(res, fname);
    return rets;
}


efp_result _efp_set_frag_coordinates(efp* efp, size_t frag_idx, efp_coord_type ctype, py::list coord) {
    enum efp_result res;

    double *ccoords = NULL;
    ccoords = new double[12];  // room for xyzabc (6), points (9), or rotmat (12)
    double *pcoords = ccoords;
    for (auto itm : coord)
        *pcoords++ = itm.cast<double>();

    res = efp_set_frag_coordinates(efp, frag_idx, ctype, ccoords);
    return res;
}


efp_result _efp_set_point_charges(efp* efp, size_t n_ptc, py::list ptc, py::list xyz) {
    enum efp_result res;

    double *cptc = NULL;
    cptc = new double[n_ptc];
    double *pptc = cptc;
    for (auto itm : ptc)
        *pptc++ = itm.cast<double>();

    double *cxyz = NULL;
    cxyz = new double[3 * n_ptc];
    double *pxyz = cxyz;
    for (auto itm : xyz)
        *pxyz++ = itm.cast<double>();

    res = efp_set_point_charges(efp, n_ptc, cptc, cxyz);
    return res;
}



py::tuple _efp_get_frag_charge(efp* efp, size_t frag_idx) {
    enum efp_result res;
    double charge = 0.0;

    res = efp_get_frag_charge(efp, frag_idx, &charge);

    py::tuple rets = py::make_tuple(res, charge);
    return rets;
}

py::tuple _efp_get_frag_multiplicity(efp* efp, size_t frag_idx) {
    enum efp_result res;
    int multiplicity = 0;

    res = efp_get_frag_multiplicity(efp, frag_idx, &multiplicity);

    py::tuple rets = py::make_tuple(res, multiplicity);
    return rets;
}

py::tuple _efp_get_gradient(efp* efp, size_t n_frag) {
    enum efp_result res;
    py::list grad;

    double *ccoords = NULL;
    ccoords = new double[6 * n_frag];
    double *pcoords = ccoords;

    res = efp_get_gradient(efp, ccoords);
    for (size_t ic = 0; ic < 6*n_frag; ++ic)
        grad.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, grad);
    return rets;
}

py::tuple _efp_get_frag_count(efp* efp) {
    enum efp_result res;
    size_t n_frag = 0;

    res = efp_get_frag_count(efp, &n_frag);

    py::tuple rets = py::make_tuple(res, n_frag);
    return rets;
}

py::tuple _efp_get_frag_atom_count(efp* efp, size_t frag_idx) {
    enum efp_result res;
    size_t n_atom = 0;

    res = efp_get_frag_atom_count(efp, frag_idx, &n_atom);

    py::tuple rets = py::make_tuple(res, n_atom);
    return rets;
}

py::tuple _efp_get_multipole_count(efp* efp) {
    enum efp_result res;
    size_t n_mult = 0;

    res = efp_get_multipole_count(efp, &n_mult);

    py::tuple rets = py::make_tuple(res, n_mult);
    return rets;
}

py::tuple _efp_get_multipole_coordinates(efp* efp, size_t n_mult) {
    enum efp_result res;
    py::list xyz;

    size_t dim = 3 * n_mult;
    double *ccoords = NULL;
    ccoords = new double[dim];
    double *pcoords = ccoords;

    res = efp_get_multipole_coordinates(efp, ccoords);
    for (size_t ic = 0; ic < dim; ++ic)
        xyz.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, xyz);
    return rets;
}

py::tuple _efp_get_multipole_values(efp* efp, size_t n_mult) {
    enum efp_result res;
    py::list mult;

    size_t dim = (1 + 3 + 6 + 10) * n_mult;
    double *ccoords = NULL;
    ccoords = new double[dim];
    double *pcoords = ccoords;

    res = efp_get_multipole_values(efp, ccoords);
    for (size_t ic = 0; ic < dim; ++ic)
        mult.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, mult);
    return rets;
}

py::tuple _efp_get_induced_dipole_count(efp* efp) {
    enum efp_result res;
    size_t n_dip = 0;

    res = efp_get_induced_dipole_count(efp, &n_dip);

    py::tuple rets = py::make_tuple(res, n_dip);
    return rets;
}

py::tuple _efp_get_induced_dipole_coordinates(efp* efp, size_t n_dip) {
    enum efp_result res;
    py::list xyz;

    size_t dim = 3 * n_dip;
    double *ccoords = NULL;
    ccoords = new double[dim];
    double *pcoords = ccoords;

    res = efp_get_induced_dipole_coordinates(efp, ccoords);
    for (size_t ic = 0; ic < dim; ++ic)
        xyz.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, xyz);
    return rets;
}

py::tuple _efp_get_induced_dipole_values(efp* efp, size_t n_dip) {
    enum efp_result res;
    py::list vals;

    size_t dim = 3 * n_dip;
    double *ccoords = NULL;
    ccoords = new double[dim];
    double *pcoords = ccoords;

    res = efp_get_induced_dipole_values(efp, ccoords);
    for (size_t ic = 0; ic < dim; ++ic)
        vals.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, vals);
    return rets;
}

py::tuple _efp_get_induced_dipole_conj_values(efp* efp, size_t n_dip) {
    enum efp_result res;
    py::list vals;

    size_t dim = 3 * n_dip;
    double *ccoords = NULL;
    ccoords = new double[dim];
    double *pcoords = ccoords;

    res = efp_get_induced_dipole_conj_values(efp, ccoords);
    for (size_t ic = 0; ic < dim; ++ic)
        vals.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, vals);
    return rets;
}

// TODO: probably filling would fail if res not good in a lot of these


py::tuple _efp_get_frag_atoms(efp* efp, size_t frag_idx, size_t frag_natom) {
    enum efp_result res;
    struct efp_atom atoms[frag_natom];
    py::list ats_info;

    if ((res = efp_get_frag_atoms(efp, frag_idx, frag_natom, atoms))) {
        py::tuple rets = py::make_tuple(res, ats_info);
        return rets;
    }

    for (size_t iat = 0; iat < frag_natom; ++iat) {
        py::dict at_info;
        at_info[py::str("Z")] = atoms[iat].znuc;
        at_info[py::str("mass")] = atoms[iat].mass;
        at_info[py::str("label")] = atoms[iat].label;
        at_info[py::str("x")] = atoms[iat].x;
        at_info[py::str("y")] = atoms[iat].y;
        at_info[py::str("z")] = atoms[iat].z;
        ats_info.append(at_info);
    }

    py::tuple rets = py::make_tuple(res, ats_info);
    return rets;
}

//py::dict extend_efp_get_atoms(efp* efp) {
//    enum efp_result res;
//    size_t frag_natom, natom=0;
//    double frag_chg;
//    int frag_mult;
//    py::list fr, frt, frcg, frmp, full_atoms;
//
//    py::dict mol_info;
//    mol_info[py::str("units")] = "Bohr";
//    mol_info[py::str("input_units_to_au")] = 1.0;
//    mol_info[py::str("fix_com")] = true;
//    mol_info[py::str("fix_orientation")] = true;
//    mol_info[py::str("fix_symmetry")] = "c1";
//
//    size_t frag_count = twrapped_efp_get_frag_count(efp);
//    for (int ifr = 0; ifr < frag_count; ++ifr) {
//        frag_natom = twrapped_efp_get_frag_atom_count(efp, ifr);
//        py::list f;
//        f.append(natom);
//        natom = natom + frag_natom;
//        f.append(natom);
//        fr.append(f);
//        frt.append("Real");
//        frcg.append(static_cast<int>(wrapped_efp_get_frag_charge(efp, ifr)));
//        frmp.append(wrapped_efp_get_frag_multiplicity(efp, ifr));
//
//        struct efp_atom atoms[frag_natom];
//        if ((res = efp_get_frag_atoms(efp, ifr, frag_natom, atoms))) {
//            std::string sres = "efp_get_frag_atoms: " + rts(res) + "\n";
//            throw libefpException(sres.c_str());
//        }
//
//        for (size_t iat = 0; iat < frag_natom; ++iat) {
//            py::dict at_init;
//            at_init[py::str("Z")] = atoms[iat].znuc;
//            at_init[py::str("mass")] = atoms[iat].mass;
//            at_init[py::str("label")] = atoms[iat].label;
//            at_init[py::str("x")] = atoms[iat].x;
//            at_init[py::str("y")] = atoms[iat].y;
//            at_init[py::str("z")] = atoms[iat].z;
//            full_atoms.append(at_init);
//        }
//    }
//
//    mol_info[py::str("fragments")] = fr;
//    mol_info[py::str("fragment_types")] = frt;
//    mol_info[py::str("fragment_charges")] = frcg;
//    mol_info[py::str("fragment_multiplicities")] = frmp;
//    mol_info[py::str("full_atoms")] = full_atoms;
//
//    return mol_info;
//}


py::tuple _efp_get_coordinates(efp* efp, size_t n_frag) {
    enum efp_result res;
    py::list xyzabc;

    double *ccoords = NULL;
    ccoords = new double[6 * n_frag];
    double *pcoords = ccoords;

    res = efp_get_coordinates(efp, pcoords);
    for (size_t ic = 0; ic < 6*n_frag; ++ic)
        xyzabc.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, xyzabc);
    return rets;
}


py::tuple _efp_get_frag_xyzabc(efp* efp, size_t frag_idx) {
    enum efp_result res;
    py::list xyzabc;

    double *ccoords = NULL;
    ccoords = new double[6];
    double *pcoords = ccoords;

    res = efp_get_frag_xyzabc(efp, frag_idx, pcoords);
    for (size_t ic = 0; ic < 6; ++ic)
        xyzabc.append(ccoords[ic]);

    py::tuple rets = py::make_tuple(res, xyzabc);
    return rets;
}

py::tuple _efp_get_wavefunction_dependent_energy(efp* efp) {
    enum efp_result res;
    double ene;

    res = efp_get_wavefunction_dependent_energy(efp, &ene);

    py::tuple rets = py::make_tuple(res, ene);
    return rets;
}

// stash the python E_field fn here
py::function field_fn_callback;
py::function get_field_fn_callback(void) { return field_fn_callback; }

efp_result cwrapped_field_fn(size_t n_pt, const double* xyz, double* field, void* user_data) {

    // prepare xyz: double* --> py::list
    py::list py_xyz;
    for (size_t ic = 0; ic < 3 * n_pt; ++ic)
        py_xyz.append(xyz[ic]);

    py::function fn = get_field_fn_callback();
    py::list py_field = fn(py_xyz);

    // harvest field: py::list --> double*
    double* field_p = field;
    for (auto itm : py_field)
        *field_p++ = itm.cast<double>();

    return EFP_RESULT_SUCCESS;
}

void _efp_set_electron_density_field_fn(efp* efp, py::function fn) {
    field_fn_callback = fn;
    efp_set_electron_density_field_fn(efp, cwrapped_field_fn);
}

py::function dummy;
void _clear_electron_density_field_fn(efp* efp) { field_fn_callback = dummy; }

void _clean(efp* efp) {
    _clear_electron_density_field_fn(efp);
    efp_shutdown(efp);
}


//      //  py::list coord --> double* ccoords
//          double* ccoords = NULL;
//          ccoords = new double[12];  // room for xyzabc (6), points (9), or rotmat (12)
//          double* pcoords = ccoords;
//          for (auto itm : coord)
//              *pcoords++ = itm.cast<double>();
//
//      //  double* coords --> py::list
//          py::list coord;
//
//          double* ccoords = NULL;
//          ccoords = new double[12];
//          double* pcoords = ccoords;
//
//          double* ccoords = new double[12];  // alt
//
//          res = efp_get_multipole_coordinates(efp, ccoords);
//          for (size_t ic = 0; ic < 12; ++ic)
//              coord.append(ccoords[ic]);



PYBIND11_MODULE(core, m) {
    m.doc() = "Python wrapping of Parallel implementation of the Effective Fragment Potential (EFP) method";

    m.attr("__version__") = py::str("1.1");
    py::exception<libefpException>(m, "libefpException");

    py::enum_<efp_result>(m, "efp_result", "Result of a libefp operation")
        .value("EFP_RESULT_SUCCESS", EFP_RESULT_SUCCESS)                      /* Operation was successful. */
        .value("EFP_RESULT_FATAL", EFP_RESULT_FATAL)                          /* Fatal error has occurred. */
        .value("EFP_RESULT_NO_MEMORY", EFP_RESULT_NO_MEMORY)                  /* Insufficient memory. */
        .value("EFP_RESULT_FILE_NOT_FOUND", EFP_RESULT_FILE_NOT_FOUND)        /* File not found. */
        .value("EFP_RESULT_SYNTAX_ERROR", EFP_RESULT_SYNTAX_ERROR)            /* Syntax error. */
        .value("EFP_RESULT_UNKNOWN_FRAGMENT", EFP_RESULT_UNKNOWN_FRAGMENT)    /* Unknown EFP fragment. */
        .value("EFP_RESULT_POL_NOT_CONVERGED", EFP_RESULT_POL_NOT_CONVERGED)  /* Polarization SCF procedure did not converge. */
        .export_values();

    py::enum_<efp_term>(m, "efp_term", py::arithmetic(), "Flags to specify EFP energy terms")
        .value("EFP_TERM_ELEC", EFP_TERM_ELEC)                                /* EFP/EFP electrostatics. */
        .value("EFP_TERM_POL", EFP_TERM_POL)                                  /* EFP/EFP polarization. */
        .value("EFP_TERM_DISP", EFP_TERM_DISP)                                /* EFP/EFP dispersion. */
        .value("EFP_TERM_XR", EFP_TERM_XR)                                    /* EFP/EFP exchange repulsion. */
        .value("EFP_TERM_CHTR", EFP_TERM_CHTR)                                /* EFP/EFP charge transfer, reserved for future. */
        .value("EFP_TERM_AI_ELEC", EFP_TERM_AI_ELEC)                          /* Ab initio/EFP electrostatics. */
        .value("EFP_TERM_AI_POL", EFP_TERM_AI_POL)                            /* Ab initio/EFP polarization. */
        .value("EFP_TERM_AI_DISP", EFP_TERM_AI_DISP)                          /* Ab initio/EFP dispersion, reserved for future. */
        .value("EFP_TERM_AI_XR", EFP_TERM_AI_XR)                              /* Ab initio/EFP exchange repulsion, reserved for future. */
        .value("EFP_TERM_AI_CHTR", EFP_TERM_AI_CHTR)                          /* Ab initio/EFP charge transfer, reserved for future. */
        .export_values();

    py::enum_<efp_disp_damp>(m, "efp_disp_damp", "Fragment-fragment dispersion damping type")
        .value("EFP_DISP_DAMP_OVERLAP", EFP_DISP_DAMP_OVERLAP)                /* Overlap-based damping (default). */
        .value("EFP_DISP_DAMP_TT", EFP_DISP_DAMP_TT)                          /* Tang-Toennies damping. */
        .value("EFP_DISP_DAMP_OFF", EFP_DISP_DAMP_OFF)                        /* No dispersion damping. */
        .export_values();

    py::enum_<efp_elec_damp>(m, "efp_elec_damp", "Fragment-fragment electrostatic damping type")
        .value("EFP_ELEC_DAMP_SCREEN", EFP_ELEC_DAMP_SCREEN)                  /* SCREEN-controlled damping (default). */
        .value("EFP_ELEC_DAMP_OVERLAP", EFP_ELEC_DAMP_OVERLAP)                /* Overlap-based damping. */
        .value("EFP_ELEC_DAMP_OFF", EFP_ELEC_DAMP_OFF)                        /* No electrostatic damping. */
        .export_values();

    py::enum_<efp_pol_damp>(m, "efp_pol_damp", "Fragment-fragment polarization damping type")
        .value("EFP_POL_DAMP_TT", EFP_POL_DAMP_TT)                            /* Tang-Toennies like damping (default). */
        .value("EFP_POL_DAMP_OFF", EFP_POL_DAMP_OFF)                          /* No polarization damping. */
        .export_values();

    py::enum_<efp_coord_type>(m, "efp_coord_type", "Describes the way fragment coordinates are specified")
        .value("EFP_COORD_TYPE_XYZABC", EFP_COORD_TYPE_XYZABC)                /* Coordinates of center of mass of a fragment and Euler angles. */
        .value("EFP_COORD_TYPE_POINTS", EFP_COORD_TYPE_POINTS)                /* Coordinates of three points belonging to a fragment. */
        .value("EFP_COORD_TYPE_ROTMAT", EFP_COORD_TYPE_ROTMAT)                /* Coordinates of fragment center of mass and its rotation matrix. */
        .export_values();


    py::enum_<efp_pol_driver>(m, "efp_pol_driver", "Driver used for solving polarization equations")
        .value("EFP_POL_DRIVER_ITERATIVE", EFP_POL_DRIVER_ITERATIVE)          /* Iterative solution of polarization equations. */
        .value("EFP_POL_DRIVER_DIRECT", EFP_POL_DRIVER_DIRECT)                /* Direct solution of polarization equations. */
        .export_values();

    py::class_<efp_opts>(m, "efp_opts", "Options controlling EFP computation")
        .def(py::init())
        .def_readwrite("terms", &efp_opts::terms)                             /* Specifies which energy terms to compute. */
        .def_readwrite("disp_damp", &efp_opts::disp_damp)                     /* Dispersion damping type (see #efp_disp_damp). */
        .def_readwrite("elec_damp", &efp_opts::elec_damp)                     /* Electrostatic damping type (see #efp_elec_damp). */
        .def_readwrite("pol_damp", &efp_opts::pol_damp)                       /* Polarization damping type (see #efp_pol_damp). */
        .def_readwrite("pol_driver", &efp_opts::pol_driver)                   /* Driver used to find polarization induced dipoles. */
        .def_readwrite("enable_pbc", &efp_opts::enable_pbc)                   /* Enable periodic boundary conditions if nonzero. */
        .def_readwrite("enable_cutoff", &efp_opts::enable_cutoff)             /* Enable frag-fraginteraction cutoff if nonzero. */
        .def_readwrite("swf_cutoff", &efp_opts::swf_cutoff);                  /* Cutoff distance for frag-frag interactions. */

    py::class_<efp_energy>(m, "efp_energy", "EFP energy terms")
        .def(py::init())
        .def_readwrite("electrostatic", &efp_energy::electrostatic)           /* EFP/EFP electrostatic energy. */
        .def_readwrite("charge_penetration", &efp_energy::charge_penetration) /* Charge penetration energy from
                                                                                 overlap-based electrostatic damping.
                                                                                 Zero if overlap-based damping is turned
                                                                                 off. */
        .def_readwrite("electrostatic_point_charges", &efp_energy::electrostatic_point_charges)
                                                                              /* Interaction energy of EFP electrostatics
                                                                                 with point charges. */
        .def_readwrite("polarization", &efp_energy::polarization)             /* All polarization energy goes here.
                                                                                 Polarization is computed self-consist-
                                                                                 ently so it can't be separated into
                                                                                 EFP/EFP and AI/EFP parts. */
        .def_readwrite("dispersion", &efp_energy::dispersion)                 /* EFP/EFP dispersion energy. */
        .def_readwrite("ai_dispersion", &efp_energy::ai_dispersion)           /* AI/EFP dispersion energy. */
        .def_readwrite("exchange_repulsion", &efp_energy::exchange_repulsion) /* EFP/EFP exchange-repulsion energy. */
        .def_readwrite("total", &efp_energy::total);                          /* Sum of all the above energy terms. */

    py::class_<efp_atom>(m, "efp_atom", "EFP atom info")
        .def(py::init())
        //.def_readonly("label", &efp_atom::label)  // char label[32];        /* Atom label. */
        .def_readwrite("x", &efp_atom::x)                                     /* X coordinate of atom position. */
        .def_readwrite("y", &efp_atom::y)                                     /* Y coordinate of atom position. */
        .def_readwrite("z", &efp_atom::z)                                     /* Z coordinate of atom position. */
        .def_readwrite("mass", &efp_atom::mass)                               /* Atom mass. */
        .def_readwrite("znuc", &efp_atom::znuc);                              /* Nuclear charge. */

    py::class_<efp, std::unique_ptr<efp, py::nodelete>>(m, "efp", "Main libefp opaque structure")
        .def(py::init(&efp_create), "Creates a new efp object via `efp_create`")
        .def("banner", _efp_banner, "Gets a human readable banner string with information about the library")
//        .def("set_error_log", &efp_set_error_log, "Sets the error log callback function")
        .def("_efp_set_opts", &efp_set_opts, "Wrapped set computation options")
        .def("_efp_get_opts", &efp_get_opts, "Gets currently set computation options")
        .def("_efp_add_potential", &efp_add_potential, "Wrapped adds EFP potential from full file path")
        .def("_efp_add_fragment", &efp_add_fragment, "Wrapped adds a new fragment to the EFP subsystem")
        .def("_efp_prepare", &efp_prepare, "Wrapped prepares the calculation")
//        .def("skip_fragments", &efp_skip_fragments, "Skip interactions between the fragments *arg0* and *arg1* inclusive if *arg2*")
        .def("set_electron_density_field_fn", _efp_set_electron_density_field_fn, "Sets the callback function which computes electric field from electrons in ab initio subsystem")
        .def("clear_electron_density_field_fn", _clear_electron_density_field_fn, "Detaches callback function from EFP instance (necessary for destruction. Called by clean or call alongside shutdown")
        .def("_efp_set_point_charges", _efp_set_point_charges, "Wrapped setup arbitrary point charges of magnitude at locations interacting with EFP subsystem")
//        .def("get_point_charge_count", &efp_get_point_charge_count, "Gets the number of currently set point charges and return it in *arg0*")
//        .def("get_point_charge_values", &efp_get_point_charge_values, "Gets values of currently set point charges and returns them in *arg1*")
//        .def("set_point_charge_values", &efp_set_point_charge_values, "Sets values of point charges *arg0*")
//        .def("get_point_charge_coordinates", &efp_get_point_charge_coordinates, "Gets coordinates of currently set point charges and returns them in *arg1*")
//        .def("set_point_charge_coordinates", &efp_set_point_charge_coordinates, "Sets coordinates *arg0* of point charges")
//        .def("get_point_charge_gradient", &efp_get_point_charge_gradient, "Gets gradient on point charges from EFP subsystem and returns them in *arg1*")
//        .def("set_coordinates", &efp_set_coordinates, "Update positions and orientations of all fragments with types in array *arg0* and returns them in *arg1*")
        .def("_efp_set_periodic_box", &efp_set_periodic_box, "Sets up periodic box size of *arg0* by *arg1* by *arg2*")
        .def("_efp_get_periodic_box", _efp_get_periodic_box, "Gets periodic box size of *arg0* by *arg1* by *arg2*")
        .def("_efp_set_frag_coordinates", _efp_set_frag_coordinates, "Wrapped updates position and orientation of the specified effective fragment with type")
        .def("_efp_get_coordinates", _efp_get_coordinates, "Wrapped gets center of mass positions and Euler angles of the effective fragments")
        .def("_efp_get_frag_xyzabc", _efp_get_frag_xyzabc, "Wrapped gets center of mass position and Euler angles on fragment")
//        .def("get_stress_tensor", &efp_get_stress_tensor, "Gets the stress tensor and returns it in *arg0*")
//        .def("get_ai_screen", &efp_get_ai_screen, "Gets the ab initio screening parameters on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("set_orbital_energies", &efp_set_orbital_energies, "Sets ab initio orbital energies to *efp0* number core orbitals, *efp1* number active orbitals, *efp2* number virtual orbitals, *efp3* array of orbital energies")
//        .def("set_dipole_integrals", &efp_set_dipole_integrals, "Sets ab initio dipole integrals to  *efp0* number core orbitals, *efp1* number active orbitals, *efp2* number virtual orbitals, *efp3* dipole integral matrices")
        .def("_efp_get_wavefunction_dependent_energy", _efp_get_wavefunction_dependent_energy, "Wrapped updates wavefunction-dependent energy terms")
        .def("_efp_compute", &efp_compute, py::arg("do_gradient") = false, "Perform the EFP computation, optionally doing gradient")
        .def("_efp_get_frag_charge", _efp_get_frag_charge, "Gets total charge on fragment")
        .def("_efp_get_frag_multiplicity", _efp_get_frag_multiplicity, "Gets spin multiplicity on fragment")
//        .def("get_frag_multipole_count", &efp_get_frag_multipole_count, "Gets number of electrostatic multipole points on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("_efp_get_multipole_count", _efp_get_multipole_count, "Gets total number of multipoles from EFP electrostatics")
        .def("_efp_get_multipole_coordinates", _efp_get_multipole_coordinates, "Gets coordinates of electrostatics multipoles")
        .def("_efp_get_multipole_values", _efp_get_multipole_values, "Gets electrostatics multipoles from EFP fragments")
        .def("_efp_get_induced_dipole_count", _efp_get_induced_dipole_count, "Gets the number of polarization induced dipoles")
        .def("_efp_get_induced_dipole_coordinates", _efp_get_induced_dipole_coordinates, "Gets coordinates of induced dipoles")
        .def("_efp_get_induced_dipole_values", _efp_get_induced_dipole_values, "Gets values of polarization induced dipoles")
        .def("_efp_get_induced_dipole_conj_values", _efp_get_induced_dipole_conj_values, "Gets values of polarization conjugated induced dipoles")
//        .def("get_lmo_count", &efp_get_lmo_count, "Gets the number of LMOs in a fragment and returns it in *arg0*")
//        .def("get_lmo_coordinates", &efp_get_lmo_coordinates, "Gets coordinates of LMO centroids on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("get_xrfit", &efp_get_xrfit, "Gets parameters of fitted exchange-repulsion on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("_efp_get_energy", &efp_get_energy, "Gets computed energy components")
        .def("_efp_get_gradient", _efp_get_gradient, "Gets computed EFP energy gradient")
//        .def("get_atomic_gradient", &efp_get_atomic_gradient, "Gets computed EFP energy gradient on individual atoms and returns it in *arg0*")
        .def("_efp_get_frag_count", _efp_get_frag_count, "Gets the number of fragments in this computation")
        .def("_efp_get_frag_name", _efp_get_frag_name, "Gets the name of the specified effective fragment")
//        .def("get_frag_mass", &efp_get_frag_mass, "Gets total mass on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("get_frag_inertia", &efp_get_frag_inertia, "Gets fragment principal moments of inertia on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("_efp_get_frag_atom_count", _efp_get_frag_atom_count, "Gets the number of atoms on fragment")
        .def("_efp_get_frag_atoms", &efp_get_frag_atoms, "Wrapped get atoms comprising the specified 0-indexed fragment")
        .def("_efp_get_frag_atoms", _efp_get_frag_atoms, "Gets atoms comprising the specified 0-indexed fragment")
//        .def("get_electric_field", &efp_get_electric_field, "Gets electric field for a point on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("torque_to_derivative", &efp_torque_to_derivative, "Convert rigid body torque *arg1* to derivatives *arg2* of energy by Euler angles *arg0*")
//        .def("result_to_string", &efp_result_to_string, "Result value to be converted to string");
        .def("clean", &_clean, "Preferred destructor combining libefp::efp_shutdown and field_fn release")
        .def("shutdown", &efp_shutdown, "Release all resources used by this EFP");
}

// Unwrapped
// * efp_print_banner
// * efp_set_electron_density_field_user_data  # what is this?


