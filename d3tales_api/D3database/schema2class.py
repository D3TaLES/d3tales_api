import python_jsonschema_objects as pjs
from urllib import request
import json


class Schema2Class:
    # Copyright 2021, University of Kentucky
    def __init__(self, database=None, schema_name=None, schema_directory=None, named_only=False):

        # fetch schema
        self.database = database
        self.schema_name = '{}/{}'.format(schema_directory, schema_name) if schema_directory else schema_name
        schema_url = "https://raw.githubusercontent.com/D3TaLES/schema/main/schema_{}/{}.schema.json".format(
            self.database, self.schema_name).replace("robotics_", "")
        # print(schema_url)
        response = request.urlopen(schema_url)
        self.schema = json.loads(response.read().decode())
        # generating classes
        builder = pjs.ObjectBuilder(self.schema)
        ns = builder.build_classes(named_only=named_only)

        # get all name space
        for name_space in dir(ns):
            setattr(self, name_space, getattr(ns, name_space))

        # # highest-level name space for validation
        # self.high_level = getattr(ns, schema_name.title().replace("_", ""))

        # required values
        if self.schema.get("required",):
            self.required = self.schema.get("required")
        else:
            self.required = None

    def md_table(self, prop_level="properties"):
        out_text = "| <div style='width:250px'>Property</div> | Description |  Data Type  |\n| :---------- | :------------- | :------------- |\n"
        props = self.schema.get(prop_level, {})
        for prop, prop_specifications in props.items():
            description = prop_specifications.get("description")
            data_type = prop_specifications.get("type")
            out_text += "| `{}` | {} | `{}` |\n".format(prop, description, data_type)
        return out_text
