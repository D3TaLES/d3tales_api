from d3tales_api.Calculators.utils import *
from d3tales_api.Calculators.periodic_table import periodictable
import cclib
import pandas as pd
import dbstep.Dbstep as db
from rdkit.Chem import rdMolAlign
from pymatgen.core.sites import Site
from pymatgen.core.structure import Molecule
from ocelot.routines.conformerparser import pmgmol_to_rdmol


class EnergyDiffCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Solvation energy calculation

        :connection points
        "energy_final" : energy final (default = eV)
        "energy_inital" : energy initial (default = eV)

        :return: solvation energy in units eV
        """
        conns = self.make_connections(data)
        energy = (unit_conversion(conns["energy_final"], default_unit='eV') - unit_conversion(conns["energy_inital"],
                                                                                              default_unit='eV'))
        return float(np.format_float_scientific(energy, precision=precision))


class ReorganizationCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Solvation energy calculation

        :connection points
        "gs_opt" : ground state optimized energy (default = eV)
        "ion_opt" : ion optimized energy (default = eV)
        "gs_energy" : ground state energy at ion geometry (default = eV)
        "ion_energy" : ion energy at ground state geometry (default = eV)

        :return: reorganization energy in units eV
        """
        conns = self.make_connections(data)
        lambda1 = (unit_conversion(conns["ion_energy"], default_unit='eV') - unit_conversion(conns["ion_opt"],
                                                                                             default_unit='eV'))
        lambda2 = (unit_conversion(conns["gs_energy"], default_unit='eV') - unit_conversion(conns["gs_opt"],
                                                                                            default_unit='eV'))
        return float(np.format_float_scientific(lambda1 + lambda2, precision=precision))


class RelaxationCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Relaxation energy calculation from high energy geometry to optimized geometry

        :connection points
        "opt_energy" : optimized energy (default = eV)
        "energy" : energy another geometry (default = eV)

        :return: relaxation energy in units eV
        """
        conns = self.make_connections(data)
        lambdax = unit_conversion(conns["energy"], default_unit='eV') - unit_conversion(conns["opt_energy"],
                                                                                        default_unit='eV')
        return float(np.format_float_scientific(lambdax, precision=precision))


class RMSECalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Root mean squared error calculator

        :connection points
        "geom_final" : geometry final (default = A)
        "geom_inital" : geometry inital (default = A)

        :return: RMSE in units A
        """
        conns = self.make_connections(data)
        geom1 = pmgmol_to_rdmol(Molecule.from_sites([Site.from_dict(sd) for sd in conns["geom_inital"]]))[0]
        geom2 = pmgmol_to_rdmol(Molecule.from_sites([Site.from_dict(sd) for sd in conns["geom_final"]]))[0]
        try:
            rmsd = rdMolAlign.GetBestRMS(geom1, geom2)
        except:
            raise ValueError("Error finding RMSE")
        return float(np.format_float_scientific(rmsd, precision=precision))


class DeltaGSolvCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Root mean squared error calculator

        :connection points
        "init_eng" : inital energy (default = eV)
        "init_corr" : initial entropy correction (default = eV)
        "init_eng_solv" : initial energy of solvation (default = eV)
        "fin_eng" : final energy (default = eV)
        "fin_corr" : final entropy correction (default = eV)
        "fin_eng_solv" : final energy of solvation (default = eV)

        :return: RMSE in units A
        """
        conns = self.make_connections(data)

        g_gas_init = unit_conversion(conns["init_eng"], default_unit='eV') + unit_conversion(conns["init_corr"],
                                                                                             default_unit='eV')
        g_gas_fin = unit_conversion(conns["fin_eng"], default_unit='eV') + unit_conversion(conns["fin_corr"],
                                                                                           default_unit='eV')

        # entropy correction cancels out because the species are the same
        delta_g_init_solv = unit_conversion(conns["init_eng_solv"], default_unit='eV') - unit_conversion(
            conns["init_eng"], default_unit='eV')
        delta_g_fin_solv = unit_conversion(conns["fin_eng_solv"], default_unit='eV') - unit_conversion(conns["fin_eng"],
                                                                                                       default_unit='eV')

        delta_g = g_gas_fin - g_gas_init + delta_g_fin_solv - delta_g_init_solv

        return float(np.format_float_scientific(delta_g, precision=precision))


class RedoxPotentialCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Root mean squared error calculator

        :connection points
        "init_eng" : initial energy (default = eV)
        "init_corr" : initial entropy correction (default = eV)
        "init_eng_solv" : initial energy of solvation (default = eV)
        "fin_eng" : final energy (default = eV)
        "fin_corr" : final entropy correction (default = eV)
        "fin_eng_solv" : final energy of solvation (default = eV)

        "num_electrons" : number of electrons (default = 1)
        "electrode" : electrode name as str or potential as float (default = standard_hydrogen_electrode)

        :return: RMSE in units A
        """
        conns = self.make_connections(data)
        delta_g = DeltaGSolvCalc(connector=self.key_pairs).calculate(data)

        potential = -delta_g / conns.get("num_electrons", 1) + get_electrode_potential(
            conns.get("electrode", "standard_hydrogen_electrode"))

        return float(np.format_float_scientific(potential, precision=precision))


class RadBuriedVolCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Radical buried volume for the atom with the atom with the highest portion
        of spin. Uses DBSTEP.

        :connection points
        "log_file" : calculation output file. Must be readable with CCLIB
        "spin_type" : type of CCLIB spin to extract (default = mulliken)

        :return: radical buried volume in units A^3
        """
        conns = self.make_connections(data)

        cmol = cclib.io.ccopen(conns["log_file"]).parse()
        spins = cmol.atomspins[conns.get("spin_type", "mulliken")]
        all_data = pd.DataFrame({"atoms": cmol.atomnos, "spin_density": spins})
        all_data["atom_sym"] = all_data.apply(lambda x: periodictable[int(x.atoms)], axis=1)
        all_data["atom_idx"] = all_data.index
        all_data["atom_idx"] = all_data.apply(lambda x: x.atom_idx + 1, axis=1)
        all_data = all_data[all_data.atom_sym != "H"]
        all_data["fractional_spin"] = all_data["spin_density"].abs() / all_data["spin_density"].abs().sum()

        self.max_cdf = all_data.loc[all_data["fractional_spin"].idxmax()]
        rad_bur_vol = float(db.dbstep(conns["log_file"], atom1=self.max_cdf["atom_idx"], volume=True).bur_vol)

        return float(np.format_float_scientific(rad_bur_vol, precision=precision))


class RadicalSpinCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Radical spin densityy. Uses DBSTEP.

        :connection points
        "log_file" : calculation output file. Must be readable with CCLIB
        "spin_type" : type of CCLIB spin to extract (default = mulliken)

        :return: radical stability score
        """
        self.make_connections(data)
        rad_bur_obj = RadBuriedVolCalc(connector=self.key_pairs)
        rad_bur_obj.calculate(data)
        max_cdf = rad_bur_obj.max_cdf

        return float(np.format_float_scientific(max_cdf["fractional_spin"], precision=precision))


class RSSCalc(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Radical stability score. Uses DBSTEP.

        :connection points
        "log_file" : calculation output file. Must be readable with CCLIB
        "spin_type" : type of CCLIB spin to extract (default = mulliken)

        :return: radical stability score
        """
        self.make_connections(data)
        rad_bur_obj = RadBuriedVolCalc(connector=self.key_pairs)
        rad_bur_vol = rad_bur_obj.calculate(data)
        max_cdf = rad_bur_obj.max_cdf
        rss = rad_bur_vol + 50 * (1 - max_cdf["fractional_spin"])

        return float(np.format_float_scientific(rss, precision=precision))
