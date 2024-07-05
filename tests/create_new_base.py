import math
import random

import pymongo
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from d3tales_api.D3database.d3database import DBconnector, db_info
from d3tales_api.D3database.back2front import Gaus2FrontCharacterization


def initialize_new_db(front_coll, old_front_coll, limit=0):
    query = {"mol_characterization.omega.0": {"$exists": True}, "mol_characterization.omega.1": {"$exists": False}}
    init_data = list(old_front_coll.find(query, {"mol_characterization.omega": 1, "mol_info": 1}).limit(limit))

    print("Inserting {} new documents...".format(len(init_data)))
    new_inti_data = [{"_id": d["_id"], "mol_info": d["mol_info"],
                      "mol_characterization": {"omega": d["mol_characterization"]["omega"][0]}} for d in init_data]
    front_coll.insert_many(new_inti_data)


def update_doc(doc):
    mol_id = doc.get("_id")
    print("Starting {}...".format(mol_id))
    omega_q = doc.get("mol_characterization", {}).get("omega")
    omega = math.floor(omega_q["value"] * 10000) / 10000
    cond = omega_q["conditions"]
    cond.update({"tuning_parameter": omega})
    g2c = Gaus2FrontCharacterization(
        _id=mol_id,
        conditions=cond,
        verbose=0,
        all_props=True
    )
    g2c.insert_all_species()
    g2c.insert_all_mol()
    print("----- Success Inserting", mol_id, " -----")
    return g2c.mol_hashes


def try_update_hashes(doc):
        used_hashes = update_doc(doc)
        hash_coll = DBconnector(db_info.get("frontend")).get_collection("used_hashes")
        try:
            hash_coll.insert_many([{"_id": h} for h in used_hashes], ordered=False)
        except pymongo.errors.BulkWriteError:
            pass

def update_ox_ids(front_coll, old_front_coll, limit=1000):
    print("Getting ox ids...")
    ox_ids = old_front_coll.find({"mol_characterization.oxidation_potential": {"$exists": True}}).distinct("_id")
    print("Getting docs to update...")
    docs_to_update = list(
        front_coll.find({"species_characterization": {"$exists": False}}, {"mol_characterization.omega": 1}).limit(
            limit))  # mol_characterization.oxidation_potential
    random.shuffle(docs_to_update)
    print(
        f"Starting multiprocessing with {multiprocessing.cpu_count()} CPUs to insert props for {len(docs_to_update)} ids")
    for d in docs_to_update:
        if d["_id"] in ox_ids:
            try_update_hashes(d)


def update_geoms(front_coll, limit=1000):
    print("Getting docs to update...")
    docs_to_update = list(
        front_coll.find({"species_characterization.groundState.geometry": {"$exists": False},
                         "species_characterization": {"$exists": True}},
                        {"mol_characterization.omega": 1}).limit(limit))
    random.shuffle(docs_to_update)
    print(
        f"Starting multiprocessing with {multiprocessing.cpu_count()} CPUs to insert props for {len(docs_to_update)} ids")
    for d in docs_to_update:
        try_update_hashes(d)



if __name__ == "__main__":
    old_coll = DBconnector(db_info.get("frontend")).get_collection("base")
    new_coll = DBconnector(db_info.get("frontend")).get_collection("base_new")

    # initialize_new_db(new_coll, old_coll)
    # update_ox_ids(new_coll, old_coll, limit=1000)
    update_geoms(new_coll, limit=1000)



    # with ThreadPoolExecutor(max_workers=8) as executor:
    #     executor.map(try_update_hashes, docs_to_update)
