from d3tales_api.D3database.back2front import CV2Front

IDS = ['a46d18cb-2f24-40fa-abae-c286e88d08ba','26e12c06-2493-47c9-b5b2-ab83e8a95ac4','88549903-7f2e-4cf9-b87a-3872b5e1ca71',
'f4b93fc0-cc82-4fdb-9f8c-ad1006c1b968','ebb5e215-d0d2-4e8e-9377-8636ffb9278f','befbd771-c80c-4de1-9dfc-88842cbbba1d',
'be3e2f05-639d-4475-84e4-ecaf1dc4f895','74cf3be4-d4ea-47fb-9368-51b04faffa90']

IDS_2P = [
'2daff38f-8942-40c2-8702-ba33ac80b465',
'0e28a47a-9aa4-453f-9ddb-9720ee340092',
'ba0529a6-d376-4cfe-8d65-b07e87fd000b',
'b06e5f8c-51fe-4d70-9871-f109d3bc110d',
'c6b8ec65-6737-4067-aa63-d28613fff584',
'5185e32b-8515-48ca-a96e-0993c0ea168e',
'13426a53-c622-40cc-9073-8545270b3bd8',
'2d0bc067-ed30-4b47-afd1-d768aefc6295'
]

if __name__ == "__main__":
    b2f = CV2Front(id_list=IDS_2P, run_processing=False, insert=False)
    print(b2f.multi_data)
