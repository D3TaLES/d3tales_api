import pytest
from d3tales_api.Processors import parser_dft

@pytest.mark.integration
@pytest.mark.gaussian
def test_gaussian_fixture_directory_exists(data_dir):
    path=data_dir / "gaussian"
    if not path.exists(): pytest.skip("Move legacy Gaussian logs to tests/data/gaussian/")
    assert any(p.suffix.lower() in {".log", ".out", ".json"} for p in path.rglob("*"))

@pytest.mark.integration
def test_dft_parser_module_exposes_parser():
    names=[n for n in dir(parser_dft) if n.lower().startswith(("process", "parse"))]
    assert names, "No public DFT parser was found"
