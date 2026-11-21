"""Database and website adapters used by Workflows.

Database-facing operations are delegated to the canonical D3database modules.
Imports stay lazy so local workflows do not require database clients or
credentials at import time.
"""

import os

D3TALES_URL = os.environ.get("D3TALES_URL", "https://d3tales.as.uky.edu")
FRONT_COLL = "base"

__all__ = [
    "DBconnector",
    "RESTAPI",
    "d3t_get_smiles",
    "d3t_get_init_structure",
    "d3t_get_ground_state",
    "d3t_get_omega",
    "d3t_get_raw_computation",
    "d3t_submit_calculation",
    "d3t_register_molecule",
]


def __getattr__(name):
    """Keep the former adapter class imports available without eager DB imports."""
    if name == "RESTAPI":
        from d3tales_api.D3database.restapi import RESTAPI
        return RESTAPI
    if name == "DBconnector":
        from d3tales_api.D3database.d3database import DBconnector
        return DBconnector
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))


def _rest_api(**kwargs):
    rest_api_class = globals().get("RESTAPI")
    if rest_api_class is None:
        from d3tales_api.D3database.restapi import RESTAPI as rest_api_class
    return rest_api_class(**kwargs)


def _rest_get(endpoint):
    return _rest_api(method="get", endpoint=endpoint, url=D3TALES_URL, return_json=True).response


def _mol_info_field(identifier, field):
    response = _rest_get("restapi/molecules/_id={}/mol_info.{}=1".format(identifier, field))
    if response:
        return response[0].get("mol_info", {}).get(field)
    return None


def d3t_get_smiles(identifier):
    """Return the frontend SMILES for an identifier, or an empty string."""
    return _mol_info_field(identifier, "smiles") or ""


def d3t_get_init_structure(identifier):
    """Return the frontend initial structure as a pymatgen Molecule."""
    from pymatgen.core.structure import Molecule

    xyz = _mol_info_field(identifier, "init_structure")
    if not xyz:
        raise LookupError("No init_structure found in D3TaLES for {}".format(identifier))
    return Molecule.from_str(xyz, "xyz")


def d3t_get_ground_state(identifier, prop="charge"):
    """Return frontend charge/spin as an integer, or None when absent."""
    value = _mol_info_field(identifier, "groundState_{}".format(prop))
    return int(value) if value is not None else None


def d3t_get_omega(identifier):
    """Return a molecule's tuned omega value, or None when absent."""
    response = _rest_get("restapi/molecules/_id={}/mol_characterization.omega=1".format(identifier))
    try:
        return response[0]["mol_characterization"]["omega"]["value"]
    except (IndexError, KeyError, TypeError):
        return None


def d3t_get_raw_computation(_hash, fields=None):
    """Return backend computation data, or None when the record is absent."""
    endpoint = "restapi/rawdata/computation/_id={}".format(_hash)
    if fields:
        endpoint += "/{}=1".format(fields)
    response = _rest_get(endpoint)
    if response:
        return response[0].get("data", {})
    return None


def d3t_submit_calculation(zip_path, identifier, calc_type):
    """Upload a zipped calculation to the D3TaLES website."""
    submission = _rest_api(
        method="post",
        endpoint="tools/upload/computation-gaussian",
        url=D3TALES_URL,
        expected_endpoint="tools/user_uploads",
        upload_file=zip_path,
        params=dict(molecule_id=identifier, calculation_type=calc_type),
    )
    if not submission.successful:
        raise RuntimeError(
            "Calculation files not successfully submitted to {} (response: {}).".format(
                submission.endpoint, submission.response.request.url
            )
        )
    return submission


def d3t_register_molecule(smiles, names=None, origin_group="Risko", public=False):
    """Register a molecule through the canonical D3database FrontDB API."""
    from d3tales_api.D3database.d3database import FrontDB
    from d3tales_api.Processors.info_from_smiles import GenerateMolInfo

    instance = GenerateMolInfo(smiles, names=names or [], origin_group=origin_group).mol_info_dict
    frontend = FrontDB(
        schema_layer="mol_info",
        instance=instance,
        smiles=smiles,
        group=origin_group,
        public=public,
        collection_name=FRONT_COLL,
    )
    return frontend.id, instance.get("groundState_charge"), instance.get("groundState_spin")
