import pytest
from d3tales_api.Processors import parser_nlp

@pytest.mark.integration
def test_nlp_fixture_directory_exists(data_dir):
    path=data_dir / "nlp"
    if not path.exists(): pytest.skip("Move legacy NLP fixtures to tests/data/nlp/")
    assert any(p.is_file() for p in path.rglob("*"))

@pytest.mark.integration
def test_nlp_parser_module_imports():
    assert parser_nlp is not None
