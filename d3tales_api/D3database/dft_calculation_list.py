def calculation_list():
    """
    Get list of D3TaLES workflow calculation names

    :return: list of calculations
    """
    return ['wtuning',
            'opt_mol',

            'opt_groundState',
            'opt_cation1',
            'opt_cation2',
            'opt_anion1',
            'opt_anion2',
            'freq_groundState',
            'freq_cation1',
            'freq_cation2',
            'freq_anion1',
            'freq_anion2',

            'solv_opt_groundState',
            'solv_opt_cation1',
            'solv_opt_cation2',
            'solv_opt_anion1',
            'solv_opt_anion2',
            'solv_energy_gsgs',
            'solv_energy_c1c1',
            'solv_energy_c2c2',
            'solv_energy_a1a1',
            'solv_energy_a2a2',
            'solv_energy_c1gs',
            'solv_energy_gsa1',
            'solv_energy_a1gs',
            'solv_energy_gsc2',
            'solv_energy_c2gs',
            'solv_energy_gsa2',
            'solv_energy_a2gs',
            'energy_gsgs',
            'energy_c1c1',
            'energy_c2c2',
            'energy_a1a1',
            'energy_a2a2',
            'energy_gsc1',
            'energy_c1gs',
            'energy_gsa1',
            'energy_a1gs',
            'energy_gsc2',
            'energy_c2gs',
            'energy_gsa2',
            'energy_a2gs',

            'tddft_groundState',
            'tddft_cation1',
            'tddft_cation2',
            'tddft_anion1',
            'tddft_anion2',
            'solv_tddft_groundState',
            'solv_tddft_cation1',
            'solv_tddft_cation2',
            'solv_tddft_anion1',
            'solv_tddft_anion2',

            'tddft']


def calc_tuple_list():
    """
    Get list of tuples for Django website dropdowns

    :return: list of tuples
    """
    data = []
    for calc in calculation_list():
        data.append((calc, calc))
    return data


if __name__ == '__main__':
    print("Calculations lists: ", calc_tuple_list())
