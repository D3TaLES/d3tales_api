import pytest

from d3tales_api.Processors.parser_uvvis import ProcessUvVis


@pytest.mark.unit
def test_process_uvvis_csv_contract(tmp_path):
    path = tmp_path / "uvvis.csv"
    path.write_text(
        "Integration Time,100 ms\n"
        "Timestamp,2024-01-01\n"
        "Instrument,Test UV\n"
        "placeholder,placeholder\n"
        "400,0.10\n"
        "500,0.20\n",
        encoding="utf-8",
    )

    data = ProcessUvVis(str(path), mol_id="M1", metadata={"solvent": "water"})

    assert data.integration_time == "100 ms"
    assert data.date_recorded == "2024-01-01"
    assert data.absorbance_data["wavelength"] == [400.0, 500.0]
    assert data.sql_data["mol_id"] == "M1"
