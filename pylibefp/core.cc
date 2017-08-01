/*
    pylibefp/core.cc: Main binding of libefp with pybind11

    Copyright (c) 2017 The Psi4 Developers

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/


#include "efp.h"
#include "private.h"

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




std::string rts(enum efp_result res) {
    std::string sres = std::string( efp_result_to_string(res) );
    return sres;
}

std::string wrapped_efp_banner(efp* efp) {
    std::string str;
    str = std::string(efp_banner());
    return str;
}


//void opts_from_dict(efp* efp, py::dict opts_init) {
efp_result opts_from_dict(efp* efp, py::dict opts_init) {
    enum efp_result res;
    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));

    for (auto item : opts_init) {
        if (std::string(py::str(item.first)) == "pol_driver") {
            if (std::string(py::str(item.second)) == "iterative")
                opts.pol_driver = EFP_POL_DRIVER_ITERATIVE;
            else if (std::string(py::str(item.second)) == "direct")
                opts.pol_driver = EFP_POL_DRIVER_DIRECT;
            else
                throw libefpException("opt_from_dict: invalid value for pol_driver");
        }

        if (std::string(py::str(item.first)) == "elec_damp") {
            if (std::string(py::str(item.second)) == "screen")
                opts.elec_damp = EFP_ELEC_DAMP_SCREEN;
            else if (std::string(py::str(item.second)) == "overlap")
                opts.elec_damp = EFP_ELEC_DAMP_OVERLAP;
            else if (std::string(py::str(item.second)) == "off")
                opts.elec_damp = EFP_ELEC_DAMP_OFF;
            else
                throw libefpException("opt_from_dict: invalid value for elec_damp");
        }

        if (std::string(py::str(item.first)) == "pol_damp") {
            if (std::string(py::str(item.second)) == "tt")
                opts.pol_damp = EFP_POL_DAMP_TT;
            else if (std::string(py::str(item.second)) == "off")
                opts.pol_damp = EFP_POL_DAMP_OFF;
            else
                throw libefpException("opt_from_dict: invalid value for pol_damp");
        }

        if (std::string(py::str(item.first)) == "disp_damp") {
            if (std::string(py::str(item.second)) == "tt")
                opts.disp_damp = EFP_DISP_DAMP_TT;
            else if (std::string(py::str(item.second)) == "overlap")
                opts.disp_damp = EFP_DISP_DAMP_OVERLAP;
            else if (std::string(py::str(item.second)) == "off")
                opts.disp_damp = EFP_DISP_DAMP_OFF;
            else
                throw libefpException("opt_from_dict: invalid value for disp_damp");
        }

        if (std::string(py::str(item.first)) == "ai_elec")
            if (opts_init["ai_elec"].cast<bool>() == true)
                opts.terms |= EFP_TERM_AI_ELEC;
        if (std::string(py::str(item.first)) == "ai_pol")
            if (opts_init["ai_pol"].cast<bool>() == true)
                opts.terms |= EFP_TERM_AI_POL;
        //opts.terms |= EFP_TERM_AI_DISP;  // may be enabled in future libefp release
        //opts.terms |= EFP_TERM_AI_XR;    // may be enabled in future libefp release
        //opts.terms |= EFP_TERM_AI_CHTR;  // may be enabled in future libefp release

        if (std::string(py::str(item.first)) == "elec")
            if (opts_init["elec"].cast<bool>() == true)
                opts.terms |= EFP_TERM_ELEC;
        if (std::string(py::str(item.first)) == "pol")
            if (opts_init["pol"].cast<bool>() == true)
                opts.terms |= EFP_TERM_POL;
        if (std::string(py::str(item.first)) == "disp")
            if (opts_init["disp"].cast<bool>() == true)
                opts.terms |= EFP_TERM_DISP;
        if (std::string(py::str(item.first)) == "xr")
            if (opts_init["xr"].cast<bool>() == true)
                opts.terms |= EFP_TERM_XR;
        //opts.terms |= EFP_TERM_CHTR;  // may be enabled in a future libefp release

        if (std::string(py::str(item.first)) == "enable_pbc")
            opts.enable_pbc = opts_init["enable_pbc"].cast<int>();
        if (std::string(py::str(item.first)) == "enable_cutoff")
            opts.enable_cutoff = opts_init["enable_cutoff"].cast<int>();
        if (std::string(py::str(item.first)) == "swf_cutoff")
            opts.swf_cutoff = opts_init["swf_cutoff"].cast<double>();
    }

    //if ((res = efp_set_opts(efp, &opts))) {
    //    std::string sres = "efp_set_opts: " + rts(res) + "\n";
    //    throw libefpException(sres.c_str());
    //}
    res = efp_set_opts(efp, &opts);
    return res;
}


py::dict opts_to_dict_old(efp* efp) {
    enum efp_result res;
    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));
    std::string ans;
    py::dict opts_init;

    efp_get_opts(efp, &opts);

    opts_init[py::str("elec")] = (opts.terms & EFP_TERM_ELEC);
    opts_init[py::str("pol")] = (opts.terms & EFP_TERM_POL);
    opts_init[py::str("disp")] = (opts.terms & EFP_TERM_DISP);
    opts_init[py::str("xr")] = (opts.terms & EFP_TERM_XR);
    //opts_init[py::str("chtr")] = (opts.terms & EFP_TERM_CHTR);

    if (opts.elec_damp == EFP_ELEC_DAMP_SCREEN)
        ans = "screen";
    else if (opts.elec_damp == EFP_ELEC_DAMP_OVERLAP)
        ans = "screen";
    else if (opts.elec_damp == EFP_ELEC_DAMP_OFF)
        ans = "off";
    opts_init[py::str("elec_damp")] = ans;

    if (opts.pol_damp == EFP_POL_DAMP_TT)
        ans = "tt";
    else if (opts.pol_damp == EFP_POL_DAMP_OFF)
        ans = "off";
    opts_init[py::str("pol_damp")] = ans;

    if (opts.disp_damp == EFP_DISP_DAMP_TT)
        ans = "tt";
    else if (opts.disp_damp == EFP_DISP_DAMP_OVERLAP)
        ans = "overlap";
    else if (opts.disp_damp == EFP_DISP_DAMP_OFF)
        ans = "off";
    opts_init[py::str("disp_damp")] = ans;

    opts_init[py::str("ai_elec")] = (opts.terms & EFP_TERM_AI_ELEC);
    opts_init[py::str("ai_pol")] = (opts.terms & EFP_TERM_AI_POL);
    //opts_init[py::str("ai_disp")] = (opts.terms & EFP_TERM_AI_DISP);
    //opts_init[py::str("ai_xr")] = (opts.terms & EFP_TERM_AI_XR);
    //opts_init[py::str("ai_chtr")] = (opts.terms & EFP_TERM_AI_CHTR);

    if (opts.pol_driver == EFP_POL_DRIVER_ITERATIVE)
        ans = "iterative";
    else if (opts.pol_driver == EFP_POL_DRIVER_DIRECT)
        ans = "direct";
    opts_init[py::str("pol_driver")] = ans;

    return opts_init;
}

py::dict opts_to_dict(efp* efp) {
    enum efp_result res;
    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));
    std::string ans;
    py::dict opts_init;

    efp_get_opts(efp, &opts);

    opts_init[py::str("elec")] = (opts.terms & EFP_TERM_ELEC) ? true : false;
    opts_init[py::str("pol")] = (opts.terms & EFP_TERM_POL) ? true : false;
    opts_init[py::str("disp")] = (opts.terms & EFP_TERM_DISP) ? true : false;
    opts_init[py::str("xr")] = (opts.terms & EFP_TERM_XR) ? true : false;
    //opts_init[py::str("chtr")] = (opts.terms & EFP_TERM_CHTR);

    if (opts.elec_damp == EFP_ELEC_DAMP_SCREEN)
        ans = "screen";
    else if (opts.elec_damp == EFP_ELEC_DAMP_OVERLAP)
        ans = "screen";
    else if (opts.elec_damp == EFP_ELEC_DAMP_OFF)
        ans = "off";
    opts_init[py::str("elec_damp")] = ans;

    if (opts.pol_damp == EFP_POL_DAMP_TT)
        ans = "tt";
    else if (opts.pol_damp == EFP_POL_DAMP_OFF)
        ans = "off";
    opts_init[py::str("pol_damp")] = ans;

    if (opts.disp_damp == EFP_DISP_DAMP_TT)
        ans = "tt";
    else if (opts.disp_damp == EFP_DISP_DAMP_OVERLAP)
        ans = "overlap";
    else if (opts.disp_damp == EFP_DISP_DAMP_OFF)
        ans = "off";
    opts_init[py::str("disp_damp")] = ans;

    opts_init[py::str("ai_elec")] = (opts.terms & EFP_TERM_AI_ELEC) ? true : false;
    opts_init[py::str("ai_pol")] = (opts.terms & EFP_TERM_AI_POL) ? true : false;
    //opts_init[py::str("ai_disp")] = (opts.terms & EFP_TERM_AI_DISP);
    //opts_init[py::str("ai_xr")] = (opts.terms & EFP_TERM_AI_XR);
    //opts_init[py::str("ai_chtr")] = (opts.terms & EFP_TERM_AI_CHTR);

    if (opts.pol_driver == EFP_POL_DRIVER_ITERATIVE)
        ans = "iterative";
    else if (opts.pol_driver == EFP_POL_DRIVER_DIRECT)
        ans = "direct";
    opts_init[py::str("pol_driver")] = ans;

    opts_init[py::str("enable_pbc")] = opts.enable_pbc;
    opts_init[py::str("enable_cutoff")] = opts.enable_cutoff;
    opts_init[py::str("swf_cutoff")] = opts.swf_cutoff;

    return opts_init;
}


//void wrapped_efp_add_fragment(efp* efp, std::string name) {
//    enum efp_result res;
//
//    if ((res = efp_add_fragment(efp, name.c_str()))) {
//        std::string sres = "wrapped_efp_add_fragment: efp_add_fragment: " + rts(res) + "\n";
//        throw libefpException(sres.c_str());
//    }
//}

efp_result wrapped_efp_set_frag_coordinates1(efp* efp, size_t frag_idx, efp_coord_type ctype, py::list coord) {
    enum efp_result res;
//    enum efp_coord_type ctype;
//
//    if (coord_type == "xyzabc")
//        ctype = EFP_COORD_TYPE_XYZABC;
//    else if (coord_type == "points")
//        ctype = EFP_COORD_TYPE_POINTS;
//    else if (coord_type == "rotmat")
//        ctype = EFP_COORD_TYPE_ROTMAT;
//    else
//        throw libefpException("efp_coord_type not found: ");

    double *ccoords = NULL;
    ccoords = new double[12];  // room for xyzabc (6), points (9), or rotmat (12)
    double *pcoords = ccoords;
    for (auto itm : coord)
        *pcoords++ = itm.cast<double>();

    res = efp_set_frag_coordinates(efp, frag_idx, ctype, ccoords);
    return res;
}


void wrapped_efp_set_frag_coordinates(efp* efp, size_t frag_idx, std::string coord_type, py::list coord) {
    enum efp_result res;
    enum efp_coord_type ctype;

    if (coord_type == "xyzabc")
        ctype = EFP_COORD_TYPE_XYZABC;
    else if (coord_type == "points")
        ctype = EFP_COORD_TYPE_POINTS;
    else if (coord_type == "rotmat")
        ctype = EFP_COORD_TYPE_ROTMAT;
    else
        throw libefpException("efp_coord_type not found: ");

    double *ccoords = NULL;
    ccoords = new double[12];  // room for xyzabc (6), points (9), or rotmat (12)
    double *pcoords = ccoords;
    for (auto itm : coord)
        *pcoords++ = itm.cast<double>();

    if ((res = efp_set_frag_coordinates(efp, frag_idx, ctype, ccoords))) {
        std::string sres = "efp_set_frag_coordinates: " + rts(res) + "\n";
        throw libefpException(sres.c_str());
    }
}


py::dict wrapped_efp_get_energy(efp* efp) {
    enum efp_result res;
    efp_energy ene;

    if ((res = efp_get_energy(efp, &ene))) {
        std::string sres = "efp_get_energy: " + rts(res) + "\n";
        throw libefpException(sres.c_str());
    }

//    if (do_grad_) {
//        SharedMatrix smgrad(new Matrix("EFP Gradient", nfrag_, 6));
//        double ** psmgrad = smgrad->pointer();
//        if ((res = efp_get_gradient(efp_, psmgrad[0])))
//            throw PsiException("EFP::compute():efp_get_gradient(): " +
//                std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//        smgrad->print_out();
//
//        outfile->Printf("  ==> EFP Gradient <==\n\n");
//
//        for (int i=0; i<nfrag_; i++) {
//            for (int j=0; j<6; j++) {
//                outfile->Printf("%14.6lf", psmgrad[i][j]);
//            }
//            outfile->Printf("\n");
//        }
//        outfile->Printf("\n");
//
//        torque_ = smgrad;
//
//        std::shared_ptr<Wavefunction> wfn = Process::environment.legacy_wavefunction();
//        wfn->set_efp_torque(smgrad);
//    }

    py::dict energies;
    energies[py::str("electrostatic")]               = ene.electrostatic;
    energies[py::str("charge_penetration")]          = ene.charge_penetration;
    energies[py::str("electrostatic_point_charges")] = ene.electrostatic_point_charges;
    energies[py::str("polarization")]                = ene.polarization;
    energies[py::str("dispersion")]                  = ene.dispersion;
    energies[py::str("exchange_repulsion")]          = ene.exchange_repulsion;
    energies[py::str("total")]                       = ene.total;

    energies[py::str("elec")]                        = ene.electrostatic +
                                                       ene.charge_penetration +
                                                       ene.electrostatic_point_charges;
    energies[py::str("xr")]                          = ene.exchange_repulsion;
    energies[py::str("pol")]                         = ene.polarization;
    energies[py::str("disp")]                        = ene.dispersion;

    return energies;
}

int wrapped_efp_get_frag_multiplicity(efp* efp, size_t frag_idx) {
    enum efp_result res;
    int multiplicity=0;

    if ((res = efp_get_frag_multiplicity(efp, frag_idx, &multiplicity))) {
        std::string sres = "efp_get_frag_multiplicity: " + rts(res) + "\n";
        throw libefpException(sres.c_str());
    }
    return multiplicity;
}

double wrapped_efp_get_frag_charge(efp* efp, size_t frag_idx) {
    enum efp_result res;
    double charge=0.0;

    if ((res = efp_get_frag_charge(efp, frag_idx, &charge))) {
        std::string sres = "efp_get_frag_charge: " + rts(res) + "\n";
        throw libefpException(sres.c_str());
    }
    return charge;
}

size_t wrapped_efp_get_frag_count(efp* efp) {
    enum efp_result res;
    size_t n=0;

    if ((res = efp_get_frag_count(efp, &n))) {
        std::string sres = "efp_get_frag_count: " + rts(res) + "\n";
        throw libefpException(sres.c_str());
    }
    return n;
}

size_t wrapped_efp_get_frag_atom_count(efp* efp, size_t frag_idx) {
    enum efp_result res;
    size_t n=0;

    if ((res = efp_get_frag_atom_count(efp, frag_idx, &n))) {
        std::string sres = "efp_get_frag_atom_count: " + rts(res) + "\n";
        throw libefpException(sres.c_str());
    }
    return n;
}

py::dict extend_efp_get_atoms(efp* efp) {
    enum efp_result res;
    size_t frag_natom, natom=0;
    double frag_chg;
    int frag_mult;
    py::list fr, frt, frcg, frmp, full_atoms;

    py::dict mol_info;
    mol_info[py::str("units")] = "Bohr";
    mol_info[py::str("input_units_to_au")] = 1.0;
    mol_info[py::str("fix_com")] = true;
    mol_info[py::str("fix_orientation")] = true;
    mol_info[py::str("fix_symmetry")] = "c1";

    size_t frag_count = wrapped_efp_get_frag_count(efp);
    for (int ifr = 0; ifr < frag_count; ++ifr) {
        frag_natom = wrapped_efp_get_frag_atom_count(efp, ifr);
        py::list f;
        f.append(natom);
        natom = natom + frag_natom;
        f.append(natom);
        fr.append(f);
        frt.append("Real");
        frcg.append(wrapped_efp_get_frag_charge(efp, ifr));
        frmp.append(wrapped_efp_get_frag_multiplicity(efp, ifr));

        struct efp_atom atoms[frag_natom];
        if ((res = efp_get_frag_atoms(efp, ifr, frag_natom, atoms))) {
            std::string sres = "efp_get_frag_atoms: " + rts(res) + "\n";
            throw libefpException(sres.c_str());
        }

        for (size_t iat = 0; iat < frag_natom; ++iat) {
            py::dict at_init;
            at_init[py::str("Z")] = atoms[iat].znuc;
            at_init[py::str("mass")] = atoms[iat].mass;
            at_init[py::str("label")] = atoms[iat].label;
            at_init[py::str("x")] = atoms[iat].x;
            at_init[py::str("y")] = atoms[iat].y;
            at_init[py::str("z")] = atoms[iat].z;
            full_atoms.append(at_init);
        }
    }

    mol_info[py::str("fragments")] = fr;
    mol_info[py::str("fragment_types")] = frt;
    mol_info[py::str("fragment_charges")] = frcg;
    mol_info[py::str("fragment_multiplicities")] = frmp;
    mol_info[py::str("full_atoms")] = full_atoms;

    return mol_info;
}


py::list wrapped_efp_get_coordinates(efp* efp) {
    enum efp_result res;

    size_t nfr = wrapped_efp_get_frag_count(efp);

    double *ccoords = NULL;
    ccoords = new double[6 * nfr];
    double *pcoords = ccoords;
//    for (auto itm : coord)
//        *pcoords++ = itm.cast<double>();

    if ((res = efp_get_coordinates(efp, pcoords))) {
        std::string sres = "efp_get_coordinates: " + rts(res) + "\n";
        throw libefpException(sres.c_str());
    }
}

py::list wrapped_efp_get_frag_xyzabc(efp* efp, size_t frag_idx) {
    enum efp_result res;
//    //std::vector<double> coords;
    py::list coords;

//    if ((res = efp_get_frag_xyzabc(efp, frag_idx, &coords[0]))) {
//        std::string sres = "efp_get_frag_xyzabc: " + rts(res) + "\n";
//        throw libefpException(sres.c_str());
//    }

    return coords;
}

std::string efp_opts_summary_old(efp* efp) {
    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));

    efp_get_opts(efp, &opts);
    //sprintf(buffer, "    units %-s\n", units_ == Angstrom ? "Angstrom" : "Bohr");

    char buffer[120];
    std::stringstream ss;

    sprintf(buffer, "  ==> EFP/EFP Setup <==\n\n");
    ss << buffer;
    sprintf(buffer, "  Number of EFP fragments: %12d\n", wrapped_efp_get_frag_count(efp));
    ss << buffer;
    sprintf(buffer, "  Electrostatics enabled?: %12d\n", (opts.terms & EFP_TERM_ELEC));
    ss << buffer;
    sprintf(buffer, "  Polarization enabled?:   %12d\n", (opts.terms & EFP_TERM_POL));
    ss << buffer;
    sprintf(buffer, "  Dispersion enabled?:     %12d\n", (opts.terms & EFP_TERM_DISP));
    ss << buffer;
    sprintf(buffer, "  Exchange enabled?:       %12d\n", (opts.terms & EFP_TERM_XR));
    ss << buffer;
    sprintf(buffer, "  Charge-Transfer enabled?:%12s\n", "undefined");
    ss << buffer;

    sprintf(buffer, "  Electrostatics damping:  %12d\n", (opts.elec_damp));
    ss << buffer;
    sprintf(buffer, "  Polarization damping:    %12d\n", (opts.pol_damp));
    ss << buffer;
    sprintf(buffer, "  Dispersion damping:      %12d\n", (opts.disp_damp));
    ss << buffer;

//    if (do_qm_) {
    sprintf(buffer, "  ==> QM/EFP Setup <==\n\n");
    ss << buffer;
    sprintf(buffer, "  Number of QM fragments:  %12d\n", -1); //, nfrag_);
    ss << buffer;
    sprintf(buffer, "  Electrostatics enabled?: %12d\n", (opts.terms & EFP_TERM_AI_ELEC));
    ss << buffer;
    sprintf(buffer, "  Polarization enabled?:   %12d\n", (opts.terms & EFP_TERM_AI_POL));
    ss << buffer;
    sprintf(buffer, "  Dispersion enabled?:     %12s\n", "undefined");
    ss << buffer;
    sprintf(buffer, "  Exchange enabled?:       %12s\n", "undefined");
    ss << buffer;
    sprintf(buffer, "  Charge-Transfer enabled?:%12s\n", "undefined");
    ss << buffer;

//    outfile->Printf("  Gradient enabled?:       %12s\n", do_grad_ ? "true" : "false");
//
//    print_efp_geometry();
//
//    if (do_qm_) {
//        outfile->Printf("  ==> QM Geometry <==\n\n");
//        molecule_->print();
//    }

    return ss.str();
}

std::string efp_opts_summary(efp* efp) {
    char buffer[120];
    std::stringstream ss;

    py::dict opts = opts_to_dict(efp);
//    //sprintf(buffer, "    units %-s\n", units_ == Angstrom ? "Angstrom" : "Bohr");
//
    sprintf(buffer, "  ==> EFP/EFP Setup <==\n\n");
    ss << buffer;
    sprintf(buffer, "  Number of EFP fragments: %12d\n", wrapped_efp_get_frag_count(efp));
    ss << buffer;
    sprintf(buffer, "  Electrostatics enabled?: %12s\n", opts["elec"].cast<bool>() ? "true" : "false");
    ss << buffer;
    sprintf(buffer, "  Polarization enabled?:   %12s\n", opts["pol"].cast<bool>() ? "true" : "false");
    ss << buffer;
    sprintf(buffer, "  Dispersion enabled?:     %12s\n", opts["disp"].cast<bool>() ? "true" : "false");
    ss << buffer;
    sprintf(buffer, "  Exchange enabled?:       %12s\n", opts["xr"].cast<bool>() ? "true" : "false");
    ss << buffer;
    sprintf(buffer, "  Charge-Transfer enabled?:%12s\n", "undefined");
    ss << buffer;

    sprintf(buffer, "  Electrostatics damping:  %12s\n", opts["elec_damp"].cast<std::string>().c_str());
    ss << buffer;
    sprintf(buffer, "  Polarization damping:    %12s\n", opts["pol_damp"].cast<std::string>().c_str());
    ss << buffer;
    sprintf(buffer, "  Dispersion damping:      %12s\n", opts["disp_damp"].cast<std::string>().c_str());
    ss << buffer;
    sprintf(buffer, "  Polarization driver:     %12s\n", opts["pol_driver"].cast<std::string>().c_str());
    ss << buffer;

////    if (do_qm_) {
//    sprintf(buffer, "  ==> QM/EFP Setup <==\n\n");
//    ss << buffer;
//    sprintf(buffer, "  Number of QM fragments:  %12d\n", -1); //, nfrag_);
//    ss << buffer;
    sprintf(buffer, "  Electrostatics enabled?: %12s\n", opts["ai_elec"].cast<bool>() ? "true" : "false");
    ss << buffer;
    sprintf(buffer, "  Polarization enabled?:   %12s\n", opts["ai_pol"].cast<bool>() ? "true" : "false");
    ss << buffer;
    sprintf(buffer, "  Dispersion enabled?:     %12s\n", "undefined");
    ss << buffer;
    sprintf(buffer, "  Exchange enabled?:       %12s\n", "undefined");
    ss << buffer;
    sprintf(buffer, "  Charge-Transfer enabled?:%12s\n", "undefined");
    ss << buffer;

////    outfile->Printf("  Gradient enabled?:       %12s\n", do_grad_ ? "true" : "false");
////
////    print_efp_geometry();
////
////    if (do_qm_) {
////        outfile->Printf("  ==> QM Geometry <==\n\n");
////        molecule_->print();
////    }

    return ss.str();
}


std::string extended_efp_geometry_str(efp* efp, double units_to_bohr=1.0) {
    char buffer[120];
    std::stringstream ss;

    sprintf(buffer, "\n");
    ss << buffer;
    sprintf(buffer, "  ==> EFP Geometry <==\n\n");
    ss << buffer;
    sprintf(buffer, "    Geometry (in %s * %12.8f):\n\n", "Bohr", units_to_bohr); //, charge = %d, multiplicity = %d:\n\n",
    ss << buffer;
    sprintf(buffer, "       Center              X                  Y                   Z       \n");
    ss << buffer;
    sprintf(buffer, "    ------------   -----------------  -----------------  -----------------\n");
    ss << buffer;

    py::dict mol_info = extend_efp_get_atoms(efp);

    for (auto atm : mol_info["full_atoms"]) {
            sprintf(buffer, "    %8s%4s   %17.12lf  %17.12lf  %17.12lf\n",
                atm["label"].cast<std::string>().c_str(), "",
                atm["x"].cast<double>() * units_to_bohr,
                atm["y"].cast<double>() * units_to_bohr,
                atm["z"].cast<double>() * units_to_bohr);
            ss << buffer;
    }
    return ss.str();
}


PYBIND11_PLUGIN(core) {
    py::module m("core", "Python wrapping of Parallel implementation of the Effective Fragment Potential (EFP) method");

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

    py::enum_<efp_term>(m, "efp_term", "Flags to specify EFP energy terms")
        .value("EFP_TERM_ELEC", EFP_TERM_ELEC)                                /* EFP/EFP electrostatics. */
        .value("EFP_TERM_POL", EFP_TERM_POL)                                  /* EFP/EFP polarization. */
        .value("EFP_TERM_DISP", EFP_TERM_DISP)                                /* EFP/EFP dispersion. */
        .value("EFP_TERM_XR", EFP_TERM_XR)                                    /* EFP/EFP exchange repulsion. */
        .value("EFP_TERM_CHTR", EFP_TERM_CHTR)                                /* EFP/EFP charge transfer, reserved for future use. */
        .value("EFP_TERM_AI_ELEC", EFP_TERM_AI_ELEC)                          /* Ab initio/EFP electrostatics. */
        .value("EFP_TERM_AI_POL", EFP_TERM_AI_POL)                            /* Ab initio/EFP polarization. */
        .value("EFP_TERM_AI_DISP", EFP_TERM_AI_DISP)                          /* Ab initio/EFP dispersion, reserved for future use. */
        .value("EFP_TERM_AI_XR", EFP_TERM_AI_XR)                              /* Ab initio/EFP exchange repulsion, reserved for future use. */
        .value("EFP_TERM_AI_CHTR", EFP_TERM_AI_CHTR)                          /* Ab initio/EFP charge transfer, reserved for future use. */
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
        .def(py::init());
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
        // double electrostatic;                /* EFP/EFP electrostatic energy. */
        // double charge_penetration;           /* Charge penetration energy from overlap-based electrostatic damping. Zero if overlap-based damping is turned off. */
        // double electrostatic_point_charges;  /* Interaction energy of EFP electrostatics with point charges. */
        // double polarization;                 /* All polarization energy goes here. Polarization is computed self-consistently so it can't be separated into EFP/EFP and AI/EFP parts. */
        // double dispersion;                   /* EFP/EFP dispersion energy. */
        // double ai_dispersion;                /* AI/EFP dispersion energy. */
        // double exchange_repulsion;           /* EFP/EFP exchange-repulsion energy. */
        // double total;                        /* Sum of all the above energy terms. */

    py::class_<efp_atom>(m, "efp_atom", "EFP atom info")
        .def(py::init());
        // char label[32];                      /* Atom label. */
        // double x;                            /* X coordinate of atom position. */
        // double y;                            /* Y coordinate of atom position. */
        // double z;                            /* Z coordinate of atom position. */
        // double mass;                         /* Atom mass. */
        // double znuc;                         /* Nuclear charge. */

    py::class_<efp>(m, "efp", "Main libefp opaque structure")
        .def(py::init())
        //.def("result", &efp_result, "Callback function which is called by libefp to obtain electric field *arg2* in the number *arg0* specified points at positions *arg1* with initialized user data *arg3*")
        .def("banner", wrapped_efp_banner, "Gets a human readable banner string with information about the library")
        .def_static("create", &efp_create, "Creates a new efp object")
//        .def("opts_default", &efp_opts_default, "Gets default values of simulation options returns in *arg0*")
//        .def("set_error_log", &efp_set_error_log, "Sets the error log callback function")
        .def("set_opts", opts_from_dict, "Set computation options to *arg0*")
        .def("get_opts", opts_to_dict, "Gets currently set computation options")
        .def("summary", efp_opts_summary, "")
        .def("add_potential", &efp_add_potential, "Adds EFP potential from full file path *arg0*")
        //.def("add_fragment", wrapped_efp_add_fragment, "Adds a new fragment *arg0* to the EFP subsystem")
        .def("add_fragment", &efp_add_fragment, "Adds a new fragment *arg0* to the EFP subsystem")
        .def("raw_prepare", &efp_prepare, "Prepares the calculation")
//        .def("skip_fragments", &efp_skip_fragments, "Skip interactions between the fragments *arg0* and *arg1* inclusive if *arg2*")
//        .def("set_electron_density_field_fn", &efp_set_electron_density_field_fn, "Sets the callback function which computes electric field from electrons in ab initio subsystem to *arg0*")
//        .def("set_electron_density_field_user_data", &efp_set_electron_density_field_user_data, "Sets user data *arg0* to be passed to ::efp_electron_density_field_fn")
//        .def("set_point_charges", &efp_set_point_charges, "Setup *arg0* arbitrary point charges of magnitude *arg1* at locations *arg2* interacting with EFP subsystem")
//        .def("get_point_charge_count", &efp_get_point_charge_count, "Gets the number of currently set point charges and return it in *arg0*")
//        .def("get_point_charge_values", &efp_get_point_charge_values, "Gets values of currently set point charges and returns them in *arg1*")
//        .def("set_point_charge_values", &efp_set_point_charge_values, "Sets values of point charges *arg0*")
//        .def("get_point_charge_coordinates", &efp_get_point_charge_coordinates, "Gets coordinates of currently set point charges and returns them in *arg1*")
//        .def("set_point_charge_coordinates", &efp_set_point_charge_coordinates, "Sets coordinates *arg0* of point charges")
//        .def("get_point_charge_gradient", &efp_get_point_charge_gradient, "Gets gradient on point charges from EFP subsystem and returns them in *arg1*")
//        .def("set_coordinates", &efp_set_coordinates, "Update positions and orientations of all fragments with types in array *arg0* and returns them in *arg1*")
        .def("set_frag_coordinates", wrapped_efp_set_frag_coordinates, "Updates position and orientation of the specified effective 0-indexed fragment *arg0* of type *arg1*")
        .def("set_frag_coordinates1", wrapped_efp_set_frag_coordinates1, "Updates position and orientation of the specified effective 0-indexed fragment *arg0* of type *arg1*")
        .def("get_coordinates", wrapped_efp_get_coordinates, "Gets center of mass positions and Euler angles of the effective fragments and returns it in *arg0*")
        .def("get_frag_xyzabc", wrapped_efp_get_frag_xyzabc, "Gets center of mass position and Euler angles on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("set_periodic_box", &efp_set_periodic_box, "Sets up periodic box size of *arg0* by *arg1* by *arg2*")
//        .def("get_stress_tensor", &efp_get_stress_tensor, "Gets the stress tensor and returns it in *arg0*")
//        .def("get_ai_screen", &efp_get_ai_screen, "Gets the ab initio screening parameters on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("set_orbital_energies", &efp_set_orbital_energies, "Sets ab initio orbital energies to *efp0* number core orbitals, *efp1* number active orbitals, *efp2* number virtual orbitals, *efp3* array of orbital energies")
//        .def("set_dipole_integrals", &efp_set_dipole_integrals, "Sets ab initio dipole integrals to  *efp0* number core orbitals, *efp1* number active orbitals, *efp2* number virtual orbitals, *efp3* dipole integral matrices")
//        .def("get_wavefunction_dependent_energy", &efp_get_wavefunction_dependent_energy, "Updates wavefunction-dependent energy terms returning in *arg0*")
        .def("raw_compute", &efp_compute, py::arg("do_gradient") = false, "Perform the EFP computation, doing gradient if *arg0*")
        .def("get_frag_charge", wrapped_efp_get_frag_charge, "Gets total charge on 0-indexed fragment *arg0*")
        .def("get_frag_multiplicity", wrapped_efp_get_frag_multiplicity, "Gets spin multiplicity on 0-indexed fragment *arg0*")
//        .def("get_frag_multipole_count", &efp_get_frag_multipole_count, "Gets number of electrostatic multipole points on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("get_multipole_count", &efp_get_multipole_count, "Gets total number of multipoles from EFP electrostatics and returns it in *arg0*")
//        .def("get_multipole_coordinates", &efp_get_multipole_coordinates, "Gets coordinates of electrostatics multipoles and returns it in *arg0*")
//        .def("get_multipole_values", &efp_get_multipole_values, "Gets electrostatics multipoles from EFP fragments and returns it in *arg0*")
//        .def("get_induced_dipole_count", &efp_get_induced_dipole_count, "Gets the number of polarization induced dipoles and returns it in *arg0*")
//        .def("get_induced_dipole_coordinates", &efp_get_induced_dipole_coordinates, "Gets coordinates of induced dipoles and returns it in *arg0*")
//        .def("get_induced_dipole_values", &efp_get_induced_dipole_values, "Gets values of polarization induced dipoles and returns it in *arg0*")
//        .def("get_induced_dipole_conj_values", &efp_get_induced_dipole_conj_values, "Gets values of polarization conjugated induced dipoles and returns it in *arg0*")
//        .def("get_lmo_count", &efp_get_lmo_count, "Gets the number of LMOs in a fragment and returns it in *arg0*")
//        .def("get_lmo_coordinates", &efp_get_lmo_coordinates, "Gets coordinates of LMO centroids on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("get_xrfit", &efp_get_xrfit, "Gets parameters of fitted exchange-repulsion on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_energy", wrapped_efp_get_energy, "Gets computed energy components")
//        .def("get_gradient", &efp_get_gradient, "Gets computed EFP energy gradient and returns it in *arg0*")
//        .def("get_atomic_gradient", &efp_get_atomic_gradient, "Gets computed EFP energy gradient on individual atoms and returns it in *arg0*")
        .def("get_frag_count", wrapped_efp_get_frag_count, "Gets the number of fragments in this computation")
//        .def("get_frag_name", &efp_get_frag_name, "Gets the name of the specified 0-indexed effective fragment *arg0* and returns it in *arg2* of length *arg1*")
//        .def("get_frag_mass", &efp_get_frag_mass, "Gets total mass on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("get_frag_inertia", &efp_get_frag_inertia, "Gets fragment principal moments of inertia on 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_frag_atom_count", wrapped_efp_get_frag_atom_count, "Gets the number of atoms on 0-indexed fragment *arg0*")
        //.def("get_frag_atoms", wrapped_efp_get_frag_atoms, "Gets atoms comprising the specified 0-indexed fragment *arg0* and returns it in *arg1*")
        .def("get_atoms", extend_efp_get_atoms, "docstring")
        .def("print_geometry", extended_efp_geometry_str, py::arg("units_to_bohr") = 1.0, "docstring")
//        .def("get_electric_field", &efp_get_electric_field, "Gets electric field for a point on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("torque_to_derivative", &efp_torque_to_derivative, "Convert rigid body torque *arg1* to derivatives *arg2* of energy by Euler angles *arg0*")
//        .def("shutdown", &efp_shutdown, "Release all resources used by this EFP");
//        .def("result_to_string", &efp_result_to_string, "Result value to be converted to string");
        ;

    return m.ptr();
}

//        .def("get_coordinates", &efp_get_coordinates, "Gets center of mass positions and Euler angles of the effective fragments and returns it in *arg0*")
//        .def("get_frag_xyzabc", &efp_get_frag_xyzabc, "Gets center of mass position and Euler angles on 0-indexed fragment *arg0* and returns it in *arg1*")
        //.def("banner", &efp_banner, "Gets a human readable banner string with information about the library")
        //.def("print_banner", &efp_print_banner, "Prints libefp banner to stdout")
        //.def("add_potential", &efp_add_potential, "Adds EFP potential from a file *arg0*")
        //.def("add_fragment", &efp_add_fragment, "Adds a new fragment *arg0* to the EFP subsystem")
        //.def("set_frag_coordinates", &efp_set_frag_coordinates, "Updates position and orientation of the specified effective 0-indexed fragment *arg0* of type *arg1* and returns it in *arg2*")
        //.def("prepare", &efp_prepare, "Prepares the calculation")
//        .def("compute", &efp_compute, "Perform the EFP computation, doing gradient if *arg0*")
//        .def("get_energy", &efp_get_energy, "Gets computed energy components and returns it in *arg0*")
//        .def("get_frag_charge", &efp_get_frag_charge, "Gets total charge on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("get_frag_multiplicity", &efp_get_frag_multiplicity, "Gets spin multiplicity on 0-indexed fragment *arg0* and returns it in *arg1*")

        //.def("set_opts", &efp_set_opts, "Set computation options to *arg0*")
        //.def("get_opts", &efp_get_opts, "Gets currently set computation options returns in *arg0*")
//        .def("get_frag_atom_count", wrapped_get_frag_atom_count, "Gets the number of atoms on 0-indexed fragment *arg0* and returns it in *arg1*")
//        .def("get_frag_count", &efp_get_frag_count, "Gets the number of fragments in this computation and returns it in *arg0*")
//        .def("get_frag_atoms", &efp_get_frag_atoms, "Gets atoms comprising the specified 0-indexed fragment *arg0* and returns it in *arg1*")

