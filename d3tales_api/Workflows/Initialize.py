import os
from rdkit.Chem import MolFromSmiles, MolToSmiles
from fireworks import FiretaskBase, explicit_serialize, FWAction
from d3tales_api.Workflows.utils import (get_data_store, get_groundState, local_identifier, get_store_dir,
                                         get_local_output_dir, write_json)
from d3tales_api.Workflows.d3tales_utils import d3t_get_smiles, d3t_register_molecule
# Backward compatibility: these classes used to live in this module, and FireWorks already stored in a
# LaunchPad reference them as d3tales_api.Workflows.Initialize.<Class>. Keep this import so they still load.
from d3tales_api.Workflows.d3tales_utils import Mail, EmailStarting, EmailFinished  # noqa: F401

INIT_GEOMETRY = "opt_mol"


def _geom_spec(mol_opt_geom):
    spec = {"mol_opt_geom": mol_opt_geom}
    if mol_opt_geom:
        spec["{}_geom".format(INIT_GEOMETRY)] = mol_opt_geom
    return spec


@explicit_serialize
class MoleculeInit(FiretaskBase):
    def run_task(self, fw_spec):
        data_store = get_data_store(self, fw_spec)
        identifier = self.get('identifier', )
        smiles = self.get('smiles', )
        name = self.get('mol_name', )
        mol_opt_geom = self.get('mol_opt_geom', None)
        name_list = [name] if name else []

        if data_store == "local":
            if not smiles:
                if identifier:
                    return FWAction(update_spec={"identifier": identifier, "data_store": data_store,
                                                 **_geom_spec(mol_opt_geom)})
                raise ValueError("MoleculeInit requires a SMILES (or an identifier) when data_store='local'.")
            clean_smiles = MolToSmiles(MolFromSmiles(smiles))
            identifier = identifier or local_identifier(clean_smiles)
            gs_charge = get_groundState(identifier, clean_smiles, prop='charge', data_store=data_store)
            gs_spin = get_groundState(identifier, clean_smiles, prop='spin', data_store=data_store)

            # Molecule-level record (local stand-in for the frontend mol_info entry)
            local_output_dir = get_local_output_dir(self, fw_spec)
            if local_output_dir:
                write_json(os.path.join(get_store_dir(identifier, local_output_dir), "mol_info.json"),
                           {"_id": identifier, "mol_info": {"smiles": clean_smiles, "names": name_list,
                                                            "groundState_charge": gs_charge,
                                                            "groundState_spin": gs_spin}})
            return FWAction(update_spec={"identifier": identifier, "smiles": clean_smiles, "gs_charge": gs_charge,
                                         "gs_spin": gs_spin, "data_store": data_store,
                                         **_geom_spec(mol_opt_geom)})

            # data_store == "d3tales"
        if identifier and d3t_get_smiles(identifier):
            return FWAction(update_spec={"identifier": identifier, **_geom_spec(mol_opt_geom)})
        group = self.get('origin_group', 'Risko')
        public = self.get('public', False)
        clean_smiles = MolToSmiles(MolFromSmiles(smiles))
        new_id, gs_charge, gs_spin = d3t_register_molecule(clean_smiles, names=name_list, origin_group=group,
                                                           public=public)
        return FWAction(update_spec={"identifier": new_id, "gs_charge": gs_charge, "gs_spin": gs_spin, **_geom_spec(mol_opt_geom)})
