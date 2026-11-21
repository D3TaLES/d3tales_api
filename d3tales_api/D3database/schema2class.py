import python_jsonschema_objects as pjs
from urllib import request
from urllib.parse import urlsplit
from functools import lru_cache
from copy import deepcopy
import json

from referencing import Registry, Resource
from referencing.jsonschema import DRAFT4


_SCHEMA_HOST = "raw.githubusercontent.com"
_SCHEMA_PREFIX = "/D3TaLES/schema/"


@lru_cache(maxsize=256)
def _fetch_schema(url):
    """Download a D3TaLES schema once per process, with a network timeout."""
    parts = urlsplit(url)
    if parts.scheme != "https" or parts.netloc != _SCHEMA_HOST or not parts.path.startswith(_SCHEMA_PREFIX):
        raise ValueError("Refusing to load schema from unexpected URL: {}".format(url))
    with request.urlopen(url, timeout=30) as response:
        return json.load(response)


def _retrieve_schema(uri):
    """Provide referenced JSON schemas to python-jsonschema-objects."""
    schema = _fetch_schema(uri)
    return Resource.from_contents(schema, default_specification=DRAFT4)


class Schema2Class:
    """
    Get D3TaLES schema from GitHub and load it to a class
    Copyright 2021, University of Kentucky
    """
    def __init__(self, database=None, schema_name=None, schema_directory=None, named_only=False, schema_version=None):
        """

        :param database: str, database name
        :param schema_name: str, schema name
        :param schema_directory: str, schema directory
        :param named_only: If true, only properties with an actual title attribute will be included in the resulting namespace
        """

        # fetch schema
        self.database = database
        self.branch = schema_version or "main"
        self.schema_name = '{}/{}'.format(schema_directory, schema_name) if schema_directory else schema_name
        schema_url = "https://raw.githubusercontent.com/D3TaLES/schema/{}/schema_{}/{}.schema.json".format(
            self.branch, self.database, self.schema_name).replace("robotics_", "")
        # print(schema_url)
        self.schema = _fetch_schema(schema_url)
        # Resolve external $refs from the same versioned schema repository.
        # python_jsonschema_objects may replace JSON Schema values with TypeProxy
        # instances while generating classes. Never pass that mutated schema to
        # jsonschema validation, or cache it for subsequent Schema2Class calls.
        builder = pjs.ObjectBuilder(deepcopy(self.schema), resolver=_retrieve_schema)
        ns = builder.build_classes(named_only=named_only)
        self.registry = Registry(retrieve=_retrieve_schema)

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
        """
        Generate Markdown table for properties in a schema level

        :param prop_level: str, schema layer from which to extract table
        :return: text for Markdown table
        """
        out_text = "| <div style='width:250px'>Property</div> | Description |  Data Type  |\n| :---------- | :------------- | :------------- |\n"
        props = self.schema.get(prop_level, {})
        for prop, prop_specifications in props.items():
            description = prop_specifications.get("description")
            data_type = prop_specifications.get("type")
            out_text += "| `{}` | {} | `{}` |\n".format(prop, description, data_type)
        return out_text
