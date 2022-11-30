from d3tales_api.Processors.back2front import Gaus2FrontCharacterization
from d3tales_api.D3database.d3database import DBconnector
from d3tales_api.database_info import db_info
from tqdm import tqdm
from sys import argv

IDS = ['502dc467c4780db94cc0c324a12c2b6b', '502dc467c4780db94cc0c324a12c2b6b', 'e230ec1287c9370451f772b8b131f914', '7a4ab0600445d41d0a0c3bed4abd5ece', 'af41192d7cd9cdb39f7a1fa14019d658', '6231f0a681fbc0ad1aeba52b4505fdf0', '82072a8cc9b047775c581e82bb01fffd']  # '9f45e5da7fc27892c06a3244a011ce88'
CALC_TYPE = 'opt_cation1'
CHECK = ''  # species_characterization.groundState.globular_volume

if __name__ == "__main__":
    front_dbc = DBconnector(db_info.get("frontend"))
    front_coll = front_dbc.get_collection("base")
    back_dbc = DBconnector(db_info.get("backend"))
    back_coll = back_dbc.get_collection("computation")
    if len(argv) > 1:
        IDS = argv[1:]
    elif not IDS:
        IDS = [str(_id) for _id in back_coll.find({"calculation_type": CALC_TYPE}).distinct('_id')]

    for i in tqdm(IDS):
        # try:
            print(i)
            data = back_coll.find_one({"_id": i})
            mol_id = data['mol_id']
            if CHECK:
                if front_coll.find({"_id": mol_id}).distinct(CHECK):
                    continue
            calculation_type = data['calculation_type']
            conditions = data['data']['conditions']
            charge = data['data']['charge']

            g2f = Gaus2FrontCharacterization(
                _id=mol_id,
                calculation_type=calculation_type,
                conditions=conditions,
                charge=charge, insert=True
            )
        # except:
        #     print("!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!")
