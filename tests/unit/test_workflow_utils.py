import pytest

from d3tales_api.Workflows.utils import runtime_from_log


@pytest.mark.unit
def test_runtime_from_log_handles_fractional_values(tmp_path):
    log = tmp_path / "job.log"
    log.write_text(
        "Job cpu time:       1.5 days 2.25 hours 30.5 minutes 15.25 seconds.\n",
        encoding="utf-8",
    )
    assert runtime_from_log(str(log)) == pytest.approx(38.763)
