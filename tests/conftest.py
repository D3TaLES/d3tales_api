from __future__ import annotations
import json
from pathlib import Path
import pytest

TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / "data"
LEGACY_DATA_DIR = TESTS_DIR / "legacy" / "raw_data"

@pytest.fixture
def db_info():
    from d3tales_api.database_info import db_info

    return db_info()(db_file_warning=False)

@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA_DIR

@pytest.fixture(scope="session")
def legacy_data_dir() -> Path:
    """Temporary bridge while old raw_data is migrated into tests/data."""
    return LEGACY_DATA_DIR

@pytest.fixture
def json_fixture():
    def load(path: Path):
        return json.loads(path.read_text(encoding="utf-8"))
    return load

@pytest.fixture
def ferricyanide_records():
    return [
        {"i_p":"0.044mA", "A":"0.0314cm^2", "v":"50mV/s", "n":1, "C":"5mM"},
        {"i_p":"0.055mA", "A":"0.0314cm^2", "v":"100mV/s", "n":1, "C":"5mM"},
        {"i_p":"0.075mA", "A":"0.0314cm^2", "v":"150mV/s", "n":1, "C":"5mM"},
        {"i_p":"0.079mA", "A":"0.0314cm^2", "v":"200mV/s", "n":1, "C":"5mM"},
        {"i_p":"0.085mA", "A":"0.0314cm^2", "v":"250mV/s", "n":1, "C":"5mM"},
        {"i_p":"0.091mA", "A":"0.0314cm^2", "v":"300mV/s", "n":1, "C":"5mM"},
    ]

@pytest.fixture
def cv_connector():
    return {"i_p":"i_p", "A":"A", "v":"v", "n":"n", "C":"C"}


@pytest.fixture
def gaussian_data(data_dir, json_fixture):
    candidates = [data_dir / "gaussian" / "gaus_ex_data.json", data_dir / "gaus_ex_data.json"]
    for path in candidates:
        if path.exists():
            return json_fixture(path)
    pytest.skip("Copy gaus_ex_data.json into tests/data/gaussian/")
