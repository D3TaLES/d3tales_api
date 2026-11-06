import importlib
import pytest

MODULES=[
 "d3tales_api.Calculators.calculators",
 "d3tales_api.Calculators.plotters",
 "d3tales_api.Processors.parser_echem",
 "d3tales_api.Processors.parser_dft",
 # "d3tales_api.Processors.parser_nlp",
]

@pytest.mark.integration
@pytest.mark.parametrize("module_name", MODULES)
def test_public_modules_import(module_name):
    importlib.import_module(module_name)
