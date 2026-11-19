import pytest

from d3tales_api.Processors import parser_echem


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_chronoamperometry_microamp_current_is_not_assigned_to_time(tmp_path, monkeypatch):
    path = tmp_path / "ca.csv"
    path.write_text(
        "2024-01-01\n"
        "File = ca.csv\n"
        "Instrument Mode = Amperometric i-t\n"
        "Init E = -0.1 (V)\n"
        "Low E = -0.1 (V)\n"
        "Step = 1\n"
        "Pulse Width = 0.1 (s)\n"
        "Time/s,Current/uA\n"
        "0,10\n"
        "1,20\n"
        "2,30\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        parser_echem.CAResistanceCalculator,
        "calculate",
        lambda self, data, **kwargs: 100.0,
    )

    ca = parser_echem.ProcessChiCA(str(path), _id="test")

    assert ca.t.tolist() == [0.0, 1.0, 2.0]
    assert ca.i == pytest.approx([1e-5, 2e-5, 3e-5])


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_ca_fixture_processing_and_resistance(data_dir):
    from d3tales_api.Processors.parser_echem import ProcessChiCA

    path = data_dir / "electrochemistry" / "ca" / "CA.txt"
    ca = ProcessChiCA(str(path), _id="test")

    assert ca.data_source == "ca"
    assert ca.x_unit == "s"
    assert ca.y_unit == "A"
    assert len(ca.t) == 20000
    assert len(ca.i) == 20000
    assert ca.measured_resistance == pytest.approx(105.11562718990892)
    assert ca.measured_conductance == pytest.approx(0.009513333333333332)
    assert ca.conductivity == {"value": pytest.approx(9.513333333333332), "unit": "mS/cm"}
    assert ca.f_slp == pytest.approx(3.451e-7)
    assert ca.r_slp == pytest.approx(-7.907e-7)


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_standard_cv_fixture_processing_and_descriptors(data_dir):
    from d3tales_api.Calculators.calculators import (
        AvgEHalfCalculator,
        CVDescriptorCalculator,
        DirtyElectrodeDetector,
    )
    from d3tales_api.Processors.parser_echem import ProcessChiCV

    path = data_dir / "electrochemistry" / "cv" / "CV.txt"
    cv = ProcessChiCV(str(path), _id="test")
    entry = cv.data_dict
    entry["num_electrodes"] = 1

    connector = {"scan_data": "data.scan_data"}
    descriptor = CVDescriptorCalculator(connector=connector)
    peaks = descriptor.peaks(entry)
    assert peaks["forward"][0][0] == pytest.approx(0.43)
    assert peaks["forward"][0][1] == pytest.approx(0.0005036)
    assert peaks["reverse"][0][0] == pytest.approx(0.23)
    assert peaks["reverse"][0][1] == pytest.approx(-0.0003463)
    assert descriptor.reversibility(entry) == ["quasi-reversible"]
    assert descriptor.e_half(entry) == [pytest.approx(0.33)]
    assert descriptor.peak_splittings(entry) == [pytest.approx(0.2)]
    assert len(descriptor.middle_sweep(entry)) == 2
    assert AvgEHalfCalculator(connector={"e": "data.e_half"}).calculate([entry]) == pytest.approx(0.33)
    assert DirtyElectrodeDetector(connector=connector).calculate(entry) is True


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_standard_cv_fixture_diffusion_real_value(data_dir):
    from numpy.exceptions import RankWarning
    from d3tales_api.Calculators.calculators import CVDiffusionCalculator
    from d3tales_api.Processors.parser_echem import ProcessChiCV

    path = data_dir / "electrochemistry" / "cv" / "CV.txt"
    cv = ProcessChiCV(
        str(path),
        _id="test",
        metadata={
            "working_electrode_surface_area": "0.0314cm^2",
            "redox_mol_concentration": "5mM",
            "temperature": "298K",
        },
    )
    entry = cv.data_dict
    entry["num_electrodes"] = 1
    connector = {
        "i_p": "data.peak_current",
        "A": "data.conditions.working_electrode_surface_area",
        "v": "data.conditions.scan_rate",
        "n": "num_electrodes",
        "C": "data.conditions.redox_mol_concentration",
    }

    with pytest.warns(RankWarning):
        result = CVDiffusionCalculator(connector=connector).calculate([entry])
    assert result == pytest.approx([0.001893, 0.0004733], rel=2e-3)


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_standard_cv_series_charge_transfer_real_value(data_dir):
    import glob
    from d3tales_api.Calculators.calculators import (
        CVChargeTransferCalculator,
        CVDiffusionCalculator,
    )
    from d3tales_api.Processors.parser_echem import ProcessChiCV

    cv_dir = data_dir / "electrochemistry" / "cv"
    entries = []
    for fn in sorted(glob.glob(str(cv_dir / "0[1-5]*.csv"))):
        cv = ProcessChiCV(
            fn,
            _id="test",
            metadata={
                "working_electrode_surface_area": "0.0314cm^2",
                "redox_mol_concentration": "5mM",
                "temperature": "298K",
            },
        )
        entry = cv.data_dict
        entry["num_electrodes"] = 1
        entries.append(entry)

    if len(entries) != 5:
        pytest.skip("Five standard CV scan files are required")

    diffusion = CVDiffusionCalculator(connector={
        "i_p": "data.peak_current",
        "A": "data.conditions.working_electrode_surface_area",
        "v": "data.conditions.scan_rate",
        "n": "num_electrodes",
        "C": "data.conditions.redox_mol_concentration",
    }).calculate(entries, sci_notation=True)
    for entry in entries:
        entry["diffusion"] = float(diffusion[1])

    charge = CVChargeTransferCalculator(connector={
        "T": "data.conditions.temperature",
        "D": "diffusion",
        "v": "data.conditions.scan_rate",
        "n": "num_electrodes",
        "X": "data.peak_splittings.0",
    }).calculate(entries, sci_notation=True)
    assert float(charge) == pytest.approx(0.003328, rel=2e-3)


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_microelectrode_cv_fixture_processing_and_calculators(data_dir):
    from d3tales_api.Calculators.calculators import (
        CVChargeTransferCalculatorMicro,
        CVDiffusionCalculatorMicro,
    )
    from d3tales_api.Processors.parser_echem import ProcessChiCV

    path = data_dir / "electrochemistry" / "cv_um" / "CVUM.txt"
    cv = ProcessChiCV(
        str(path),
        _id="test",
        micro_electrodes=True,
        metadata={
            "working_electrode_radius": "0.0011/2cm",
            "redox_mol_concentration": "0.01M",
            "temperature": "293K",
            "e_ref": "0.3V",
        },
    )
    entry = cv.data_dict
    entry["num_electrodes"] = 1

    assert cv.data_source == "cv"
    assert cv.num_scans == 2
    assert len(cv.potential) == 1600
    assert cv.i_ss["value"] == pytest.approx(3.6571024390243906e-8)
    assert cv.i_ss["unit"] == "A"
    assert cv.e_half["value"] == pytest.approx(0.387)
    assert cv.e_half["unit"] == "V"

    diffusion = CVDiffusionCalculatorMicro(connector={
        "i_ss": "data.i_ss",
        "r": "data.conditions.working_electrode_radius",
        "C": "data.conditions.redox_mol_concentration",
        "n": "num_electrodes",
    }).calculate(entry)
    assert diffusion == pytest.approx(1.7228725385354527e-5)

    entry["diffusion"] = diffusion
    charge = CVChargeTransferCalculatorMicro(connector={
        "e_half": "data.e_half",
        "e_ref": "data.e_ref",
        "n": "num_electrodes",
        "T": "data.conditions.temperature",
        "D": "diffusion",
        "r": "data.conditions.working_electrode_radius",
    }).calculate(entry, sci_notation=True)
    assert float(charge) == pytest.approx(-0.0020627818102318565)


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_ircomp_fixture_processing(data_dir):
    from d3tales_api.Processors.parser_echem import ProcessChiESI

    path = data_dir / "electrochemistry" / "ircomp" / "iRComp.txt"
    ir = ProcessChiESI(str(path), _id="test")

    assert ir.data_source == "esi"
    assert len(ir.freq) == 12
    assert ir.freq[0] == pytest.approx(82520.0)
    assert ir.z_real[0] == pytest.approx(-44830.0)
    assert ir.resistance == pytest.approx(-44830.0)
    assert ir.high_f == {"value": 100000.0, "unit": "Hz"}
    assert ir.low_f == {"value": 10000.0, "unit": "Hz"}
