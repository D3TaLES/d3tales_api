from pathlib import Path
import pytest
from d3tales_api.Processors import parser_echem

@pytest.mark.integration
@pytest.mark.electrochemistry
def test_cv_fixture_directory_exists(data_dir):
    path=data_dir / "electrochemistry" / "cv"
    if not path.exists(): pytest.skip("Move legacy CV raw files to tests/data/electrochemistry/cv/")
    assert any(p.is_file() for p in path.rglob("*"))

@pytest.mark.integration
def test_echem_parser_module_exposes_parser():
    names=[n for n in dir(parser_echem) if n.lower().startswith(("process", "parse"))]
    assert names, "No public electrochemistry parser was found"
