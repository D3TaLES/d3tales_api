from d3tales_api.Processors.back2front import CV2Front

IDS = ['a46d18cb-2f24-40fa-abae-c286e88d08ba','26e12c06-2493-47c9-b5b2-ab83e8a95ac4','88549903-7f2e-4cf9-b87a-3872b5e1ca71',
'f4b93fc0-cc82-4fdb-9f8c-ad1006c1b968','ebb5e215-d0d2-4e8e-9377-8636ffb9278f','befbd771-c80c-4de1-9dfc-88842cbbba1d',
'be3e2f05-639d-4475-84e4-ecaf1dc4f895','74cf3be4-d4ea-47fb-9368-51b04faffa90']

if __name__ == "__main__":
    b2f = CV2Front(id_list=IDS, e_half_scan_rate=0.1, run_anodic=False, run_processing=True, insert=False, verbose=1)
    print(b2f.meta_dict)