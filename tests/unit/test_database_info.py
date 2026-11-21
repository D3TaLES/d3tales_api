import json
import pytest

from d3tales_api.database_info import db_info, source_groups


@pytest.mark.unit
def test_db_info_generator_reads_json_file(tmp_path):
    expected = {
        "d3tales": {"host": "mongodb://example/d3tales", "database": "d3tales"},
    }
    path = tmp_path / "db_info.json"
    path.write_text(json.dumps(expected), encoding="utf-8")
    assert db_info(str(path)) == expected


@pytest.mark.unit
def test_db_info_generator_uses_environment_variable(tmp_path, monkeypatch):
    expected = {"d3tales": {"host": "x", "database": "d3tales"}}
    path = tmp_path / "db_info.json"
    path.write_text(json.dumps(expected), encoding="utf-8")
    monkeypatch.setenv("DB_INFO_FILE", str(path))
    assert db_info() == expected


@pytest.mark.unit
def test_db_info_generator_returns_defaults_without_file(monkeypatch):
    monkeypatch.delenv("DB_INFO_FILE", raising=False)
    result = db_info(db_file_warning=False)
    assert result["d3tales"]["database"] == "d3tales"
    assert "frontend" not in result
    assert "backend" not in result


@pytest.mark.unit
def test_db_info_generator_warns_without_file(monkeypatch):
    monkeypatch.delenv("DB_INFO_FILE", raising=False)
    with pytest.warns(UserWarning, match="DB_INFO_FILE"):
        db_info(db_file_warning=True)


@pytest.mark.unit
def test_source_groups_generator_reads_json_file(tmp_path):
    expected = {"Risko": "01", "Unknown": "00"}
    path = tmp_path / "groups.json"
    path.write_text(json.dumps(expected), encoding="utf-8")
    assert source_groups(str(path)) == expected
