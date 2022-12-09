import warnings
from rdkit.Chem import MolFromSmiles, MolToSmiles
from django.forms.models import model_to_dict
from django.core.serializers.json import DjangoJSONEncoder
from d3tales_api.Calculators.calculators import *
from d3tales_api.Calculators.utils import dict2obj
from d3tales_api.D3database.d3database import FrontDB
from d3tales_api.D3database.info_from_smiles import GenerateMolInfo


def get_id(o):
    """
    Get id (either `_id` or `uuid`) from a data object
    """
    return getattr(o, '_id', None) or getattr(o, 'uuid')


class ProcessExpFlowObj:
    """
    Base class for Processing ExpFlow objects for data analysis
    """

    def __init__(self, expflow_obj, source_group, **kwargs):
        """
        :param expflow_obj: obj or dict, ExpFlow object to process
        :param source_group: str, source group for ExpFlow object
        :param kwargs: object for processing ExpFlow objects
        """
        expflow_dict = expflow_obj if isinstance(expflow_obj, dict) else json.loads(
            json.dumps(model_to_dict(expflow_obj), cls=DjangoJSONEncoder))
        self.expflow_obj = dict2obj(expflow_dict)
        self.source_group = source_group
        self.object_id = get_id(self.expflow_obj)
        self.workflow = self.expflow_obj.workflow or []
        self.reagents = self.expflow_obj.reagents or []
        self.apparatus = self.expflow_obj.apparatus or []
        self.instruments = self.expflow_obj.instruments or []
        self.kwargs = dict(**kwargs)
        self.concentration_smiles = None
        self.concentration_volume = None
        self.concentration_mass = None

    @property
    def redox_mol(self):
        """Redox molecule instance"""
        reagent_instance = [e for e in self.reagents if e.type == "redox_molecule"]
        if len(reagent_instance) == 1:
            return reagent_instance[0]
        else:
            warnings.warn(
                "Error. ExpFlow object {} has {} redox moleucle entries".format(self.object_id, len(reagent_instance)))

    @property
    def molecule_id(self):
        """Molecule ID"""
        rdkmol = MolFromSmiles(self.redox_mol.smiles)
        clean_smiles = MolToSmiles(rdkmol)
        check_id = FrontDB(smiles=clean_smiles).check_if_in_db()
        if check_id:
            return check_id
        instance = GenerateMolInfo(clean_smiles, database="frontend").mol_info_dict
        db_insertion = FrontDB(schema_layer='mol_info', instance=instance, smiles=clean_smiles,
                               group=self.source_group)
        return db_insertion.id

    @property
    def reagents_by_id(self):
        """Dictionary of reagent ids and reagent names"""
        reagents_dict = {}
        for reagent in self.reagents:
            reagent_vars = vars(reagent)
            reagent_vars.pop(get_id(reagent))
            reagents_dict[get_id(reagent)] = reagent_vars
        return reagents_dict

    @property
    def apparatus_by_id(self):
        """Dictionary of apparatus ids and apparatus names"""
        apparatus_dict = {}
        for apparatus in self.apparatus:
            apparatus_vars = vars(apparatus)
            apparatus_vars.pop(get_id(apparatus))
            apparatus_dict[get_id(apparatus)] = apparatus_vars
        return apparatus_dict

    # Get instrument name
    @property
    def instrument_name(self, instrument_idx=0):
        """Instrument name"""
        if self.instruments:
            return self.instruments[instrument_idx].name

    def property_by_action(self, action_name, parameter_idx=0):
        """
        Get molecule_id working electrode surface area
        NOTE: this only works if there is only one action of the given name

        :param action_name: str, name of action to get
        :param parameter_idx: index of parameter to get if more than one
        :return: return property as data dict
        """

        action_instance = [e for e in self.workflow if e.name == action_name]
        if len(action_instance) == 1:
            action_value = action_instance[0].parameters[parameter_idx].value
            action_unit = action_instance[0].parameters[parameter_idx].unit
            return {"value": float(action_value), "unit": action_unit}
        else:
            warnings.warn(
                "Error. ExpFlow object {} has {} instances of {}. ExpFlow submissions must have 1 instance of {}".format(
                    self.object_id, len(action_instance), action_name, action_name))

    def get_apparatus(self, apparatus_type, get_uuid=False):
        """
        Get apparatus by type
        NOTE: this only works if there is only one apparatus of the given type

        :param apparatus_type: str, type of apparatus to get
        :param get_uuid: bool, return uuid of apparatus, else return apparatus name
        :return: apparatus name
        """
        beaker_instances = [e for e in self.apparatus if e.type == apparatus_type]
        if len(beaker_instances) > 1:
            warnings.warn(
                "Error. ExpFlow object {} has {} beakers entries".format(self.object_id, len(beaker_instances)))
        elif len(beaker_instances) == 1:
            if get_uuid:
                return get_id(beaker_instances[0])
            return beaker_instances[0].name

    def get_reagents(self, reagent_type, get_uuids=False):
        """
        Get all reagents of a given type

        :param reagent_type: str, type of reagent to get
        :param get_uuids: bool, return uuid of electrode, else return electrode name
        :return: list of reagents (as dict, EX: [{"name": water, "purity": 0.95], or as uuids)
        """
        reagent_instance = [e for e in self.reagents if e.type == reagent_type]
        instance_list = []
        for instance in reagent_instance:
            if get_uuids:
                instance_list.append(get_id(instance))
            else:
                instance_list.append(dict(name=instance.name, purity=instance.purity))
        return instance_list

    @staticmethod
    def get_param_quantity(param_list):
        """
        Get parameters with mass and/or volume quantities for all parameters in a list of parameters

        :param param_list: list of parameters
        :return: list of parameters with mass and/or volume quantities
        """
        quantity_param = [q for q in param_list if q.description in ["mass", "volume"]]
        if len(quantity_param) != 1:
            warnings.warn("Error. Parameter list has {} quantity parameters".format(len(quantity_param)))
        else:
            return quantity_param[0]

    def get_concentration(self, solvent_uuids, redox_uuid, beaker_uuid):
        """
        Set temporary properties for concentration calculator

        :param solvent_uuids: list, list of solvent uuids
        :param redox_uuid: str, uuid for redox-active molecule
        :param beaker_uuid: str, uuid for beaker
        :return: concentration data dict
        """
        self.concentration_smiles = [r for r in self.reagents if get_id(r) == redox_uuid][0].smiles
        redox_actions = [a for a in self.workflow if
                         a.start_uuid == redox_uuid and a.end_uuid == beaker_uuid and a.name == 'transfer_solid']
        if len(redox_actions) != 1:
            print(
                "Error. ExpFlow object {} has {} transfer_solid actions with the same starting and ending point".format(
                    self.object_id, len(redox_actions)))
            return None
        else:
            self.concentration_mass = self.get_param_quantity(redox_actions[0].parameters)
        solvent_actions = [a for a in self.workflow if
                           a.start_uuid in solvent_uuids and a.end_uuid == beaker_uuid and a.name == 'transfer_liquid']
        if len(solvent_actions) != 1:
            print("Error. ExpFlow object {} has {} transfer_liquid actions".format(self.object_id, len(redox_actions)))
            return None
        else:
            self.concentration_volume = self.get_param_quantity(solvent_actions[0].parameters)

        # Use D3TaLES concentration calculator
        connector = {
            "smiles": "concentration_smiles",
            "weight": "concentration_mass",
            "volume": "concentration_volume"
        }
        return {"value": ConcentrationCalculator(connector=connector).calculate(self), "unit": "M"}

    @property
    def redox_mol_concentration(self):
        """Concentration of redox-active molecule"""
        return self.get_concentration(self.get_reagents("solvent", get_uuids=True), get_id(self.redox_mol),
                                      self.get_apparatus("beaker", get_uuid=True))


class ProcessExperimentRun(ProcessExpFlowObj):
    """
    Class for Processing Experiment Run object for data analysis
    """

    def get_electrode(self, electrode, get_uuid=False):
        """
        Get electrode from experiment data
        
        :param electrode: str,which electrode (working, counter, or reference)
        :param get_uuid: bool, return uuid of electrode, else return electrode name
        :return: electrode (uuid or name)
        """
        electrode_instance = [e for e in self.apparatus if e.type == "electrode_{}".format(electrode)]
        if len(electrode_instance) > 1:
            warnings.warn("Error. ExpFlow object {} has {} {} electrode entries".format(self.object_id, len(electrode_instance), electrode))
        elif len(electrode_instance) == 1:
            if get_uuid:
                return get_id(electrode_instance[0])
            return electrode_instance[0].name

    @property
    def cv_metadata(self):
        """Dictionary of CV metadata (in accordance with D3TaLES schema)"""
        metadata = dict(
            experiment_run_id=self.object_id,
            molecule_id=self.molecule_id,
            working_electrode_surface_area=self.property_by_action("measure_working_electrode_area"),
            temperature=self.property_by_action("heat") or self.property_by_action("heat_stir"),
            instrument=self.instrument_name,
            electrode_counter=self.get_electrode("counter"),
            electrode_reference=self.get_electrode("reference"),
            electrode_working=self.get_electrode("working"),
            solvent=self.get_reagents("solvent"),
            electrolyte=self.get_reagents("electrolyte"),
            ionic_liquid=self.get_reagents("ionic_liquid"),
            redox_mol_concentration=self.redox_mol_concentration,
        )
        return {k: v for k, v in metadata.items() if v}
