from d3tales_api.Calculators.calculators import *

test_dict = [{
    "file_name": "c\\users\\susan odom\\desktop\\hussein\\2021-06-22\\csv files\\13- scan rate both waves - 25 mv-s",
    "header": "",
    "note": "",
    "date_recorded": "2021-06-22T15:50:38",
    "data_points_per_scan": 1600,
    "initial_p_n": "P",
    "segment": 1,
    "sample_interval": {
        "value": 0.001,
        "unit": "V"
    },
    "quiet_time": {
        "value": 2.0,
        "unit": "sec"
    },
    "sensitivity": {
        "value": 0.0001,
        "unit": "A/V"
    },
    "comp_r": {
        "value": 1317.2,
        "unit": "ohm"
    },
    "scan_data": [],
    "forward": [
        [
            0.78,
            7.947e-05
        ],
        [
            1.027,
            0.000109
        ]
    ],
    "reverse": [
        [
            0.921,
            -9.165e-06
        ],
        [
            0.672,
            -5.184e-05
        ]
    ],
    "conditions": {
        "scan_rate": {
            "value": 0.025,
            "unit": "V/s"
        },
        "num_scans": 2,
        "initial_potential": {
            "value": 0.0,
            "unit": "V"
        },
        "high_e": {
            "value": 1.6,
            "unit": "V"
        },
        "low_e": {
            "value": 0.0,
            "unit": "V"
        },
        "data_source": "cv",
        "working_electrode": "Polycrystalline platinum electrode",
        "counter_electrode": "Platinum wire electrode",
        "reference_electrode": "Silver/Silver electrode",
        "solvent": [
            {
                "name": "Acetonitrile",
                "purity": ""
            }
        ],
        "electrolyte": [
            {
                "name": "TBAP",
                "purity": ""
            }
        ],
        "instrument": "electrochemistry__chi_660d,_chi_1100b,_pine_wavenow",
        "working_electrode_surface_area": {
            "value": 4.0,
            "unit": "cm^2"
        },
        "redox_mol_concentration": {
            "value": 0.1,
            "unit": "mol/L"
        },
        "experiment_run_id": "d3734f7c-fdd8-4392-bc50-9f926e5c4966"
    },
    "plot_data": [
        {}
    ],
    "reversibility": [
        "quasi-reversible",
        "irreversible"
    ],
    "e_half": [
        0.8505,
        0.8495
    ],
    "peak_splittings.py": [
        0.141,
        0.355
    ],
    "middle_sweep": []
}]


def get_diffusion(processed_data, num_electrons=1, anodic=True):
    for data_dict in processed_data:
        if anodic:
            data_dict["current"] = data_dict.get("forward", {})[num_electrons][1]
        else:
            data_dict["current"] = data_dict.get("reverse", {})[num_electrons][1]
        data_dict["n"] = num_electrons
    connector = {
        "i_p": "current",
        "A": "conditions.working_electrode_surface_area.value",
        "v": "conditions.scan_rate.value",
        "C": "conditions.redox_mol_concentration.value",
        "n": "n"
    }
    diffusion_cal = CVDiffusionCalculator(connector=connector)
    return diffusion_cal.calculate(processed_data)


x = get_diffusion(test_dict)
print(x)
