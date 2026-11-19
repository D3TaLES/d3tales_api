from pathlib import Path
import pytest

from d3tales_api.Processors.parser_dft import ProcessGausLog


@pytest.mark.unit
def test_gaussian_runtime_decimal_and_missing(tmp_path):
    log = tmp_path / "runtime.log"
    log.write_text("Job cpu time: 1.0 days 2.5 hours 30.0 minutes 15.25 seconds.\n", encoding="utf-8")
    parser = ProcessGausLog.__new__(ProcessGausLog)
    parser.log_path = str(log)
    assert parser.runtime == pytest.approx(27.004)


@pytest.mark.unit
def test_gaussian_runtime_missing_job_line(tmp_path):
    log = tmp_path / "runtime.log"
    log.write_text("Normal termination\n", encoding="utf-8")
    parser = ProcessGausLog.__new__(ProcessGausLog)
    parser.log_path = str(log)
    assert parser.runtime == 0


@pytest.mark.unit
def test_gaussian_tddft_excitations_parse_spin_states(tmp_path):
    log = Path(__file__).resolve().parents[2] / "data" / "gaussian" / "05XICU" / "logs" / "05XICU_tddft_cation2_HW8.log"
    result = ProcessGausLog.get_tddft_excitations(str(log))
    assert result["Singlet"]
    assert result["Triplet"]
