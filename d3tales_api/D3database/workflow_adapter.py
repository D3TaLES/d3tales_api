"""Database and website adapters used by Workflows.

Database-facing operations are delegated to the canonical D3database modules.
Imports stay lazy so local workflows do not require database clients or
credentials at import time.
"""

import os
import zipfile

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
    "d3t_direct_get_smiles",
    "d3t_direct_get_init_structure",
    "d3t_direct_get_ground_state",
    "d3t_direct_get_omega",
    "d3t_direct_get_raw_computation",
    "d3t_direct_register_molecule",
    "d3t_direct_store_calculation",
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
    """Deprecated database-backed molecule registration for legacy callers."""
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


def _direct_collection(name):
    from d3tales_api.D3database.d3database import D3Database

    return D3Database(database="d3tales", collection_name=name, validate_schema=False).coll


def d3t_direct_get_molecule(identifier):
    """Return one molecule document from the unified direct-access database."""
    return _direct_collection("base").find_one({"_id": identifier})


def d3t_direct_get_smiles(identifier):
    document = d3t_direct_get_molecule(identifier) or {}
    return document.get("mol_info", {}).get("smiles", "")


def d3t_direct_get_init_structure(identifier):
    """Return a molecule's initial XYZ string or serialized structure."""
    document = d3t_direct_get_molecule(identifier) or {}
    structure = document.get("mol_info", {}).get("init_structure")
    if not structure:
        raise LookupError("No init_structure found in D3TaLES for {}".format(identifier))
    from pymatgen.core.structure import Molecule

    if isinstance(structure, str):
        return Molecule.from_str(structure, "xyz")
    return Molecule.from_sites(structure.get("sites", structure))


def d3t_direct_get_ground_state(identifier, prop="charge"):
    document = d3t_direct_get_molecule(identifier) or {}
    value = document.get("mol_info", {}).get("groundState_{}".format(prop))
    return int(value) if value is not None else None


def d3t_direct_get_omega(identifier, engine=None):
    document = d3t_direct_get_molecule(identifier) or {}
    omega = document.get("mol_characterization", {}).get("omega", {})
    if engine is None and isinstance(omega, dict) and omega.get("value") is not None:
        return omega["value"]
    query = {"mol_id": identifier, "data.omega": {"$exists": True}}
    if engine:
        query["data.engine"] = engine
    record = _direct_collection("computation").find_one(
        query,
        {"data.omega": 1},
        sort=[("last_updated", -1)],
    )
    if record:
        return record.get("data", {}).get("omega")
    return omega.get("value") if engine is None and isinstance(omega, dict) else None


def d3t_direct_get_raw_computation(_hash, fields=None):
    projection = {fields: 1} if fields else None
    document = _direct_collection("computation").find_one({"_id": _hash}, projection)
    if not document:
        return None
    return document.get("data", {})


def d3t_direct_register_molecule(smiles, names=None, origin_group="Risko", public=False):
    """Create or find a molecule document without website or legacy database adapters."""
    from nanoid import generate
    from d3tales_api.database_info import source_groups
    from d3tales_api.Processors.info_from_smiles import GenerateMolInfo

    collection = _direct_collection("base")
    existing = collection.find_one({"mol_info.smiles": smiles})
    if existing:
        info = existing.get("mol_info", {})
        return existing["_id"], info.get("groundState_charge"), info.get("groundState_spin")

    info = GenerateMolInfo(smiles, names=names or [], origin_group=origin_group).mol_info_dict
    group_number = source_groups().get(origin_group.title(), "10")
    identifier = group_number + generate("ABCDEFGHIJKLMNOPQRSTUVWXYZ", size=4)
    while collection.count_documents({"_id": identifier}):
        identifier = group_number + generate("ABCDEFGHIJKLMNOPQRSTUVWXYZ", size=6)
    collection.insert_one({"_id": identifier, "mol_info": info, "public": bool(public)})
    return identifier, info.get("groundState_charge"), info.get("groundState_spin")


def d3t_direct_store_calculation(record, files=None, archive_path=None):
    """Insert a processed computation document and optionally archive worker files locally."""
    if not isinstance(record, dict) or not record.get("_id"):
        raise ValueError("A direct calculation record must be a mapping with an _id.")

    from d3tales_api.D3database.d3database import D3Database
    from datetime import datetime

    payload = dict(record)
    payload.pop("_id")
    payload.setdefault("last_updated", datetime.now().isoformat())
    database = D3Database(
        database="d3tales",
        collection_name="computation",
        instance=payload,
        validate_schema=False,
    )
    database.insert(record["_id"])

    existing_files = [path for path in (files or []) if os.path.isfile(path)]
    if not existing_files or not archive_path:
        return None
    os.makedirs(os.path.dirname(os.path.abspath(archive_path)), exist_ok=True)
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in existing_files:
            archive.write(path, arcname=os.path.basename(path))
    return archive_path
