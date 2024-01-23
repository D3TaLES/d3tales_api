import os
import smtplib, ssl
from d3tales_api.D3database.d3database import FrontDB
from rdkit.Chem import MolFromSmiles, MolToSmiles
from d3tales_api.D3database.info_from_smiles import GenerateMolInfo
from fireworks import FiretaskBase, explicit_serialize, FWAction


@explicit_serialize
class MoleculeInit(FiretaskBase):
    # Copyright 2021, University of Kentucky
    def run_task(self, fw_spec):
        identifier = self.get('identifier', )
        if identifier:
            return FWAction(update_spec={"identifier": identifier})

        smiles = self.get('smiles', )
        group = self.get('origin_group', 'Risko')
        name = self.get('mol_name', )
        public = self.get('public', False)

        name_list = [name] if name else []
        rdkmol = MolFromSmiles(smiles)
        clean_smiles = MolToSmiles(rdkmol)
        instance = GenerateMolInfo(clean_smiles, names=name_list, origin_group=group, database="frontend").mol_info_dict
        gs_charge = instance.get('groundState_charge')
        gs_spin = instance.get('groundState_spin')
        db_insertion = FrontDB(schema_layer='mol_info', instance=instance, smiles=clean_smiles, group=group, public=public)
        return FWAction(update_spec={"identifier": db_insertion.id, "gs_charge": gs_charge, "gs_spin": gs_spin})


class Mail:

    def __init__(self):
        self.port = 465
        self.smtp_server_domain_name = "smtp.gmail.com"
        self.sender_mail = "d3tales@gmail.com" or os.environ["EMAIL_USER"]
        self.password = "hvcfylnzlyroxniw" or os.environ["EMAIL_PASS"]

    def send(self, emails, subject, content):
        ssl_context = ssl.create_default_context()
        service = smtplib.SMTP_SSL(self.smtp_server_domain_name, self.port, context=ssl_context)
        service.login(self.sender_mail, self.password)

        for email in emails:
            result = service.sendmail(self.sender_mail, email, f"Subject: {subject}\n\n{content}")

        service.quit()


@explicit_serialize
class EmailStarting(FiretaskBase):
    # Copyright 2021, University of Kentucky
    def run_task(self, fw_spec):
        identifier = fw_spec.get("identifier", "")
        email_address = self.get('email', "")
        username = self.get('username') or "D3TaLES User"
        content = str(username)+""",
        Your submitted molecule, """+str(identifier)+""", has been approved. Calculations for this molecule have begun. You will receive another email when calculations have finished.
        
        Thank you, 
        D3TaLES Computational Team"""

        mail = Mail()
        mail.send([email_address], "D3TaLES Molecule Calculations", content)
        return FWAction(update_spec={"identifier": identifier})


@explicit_serialize
class EmailFinished(FiretaskBase):
    # Copyright 2021, University of Kentucky
    def run_task(self, fw_spec):
        identifier = fw_spec.get("identifier", "")
        email_address = self.get('email', "")
        username = self.get('username', ) or "D3TaLES User"
        content = str(username)+""",
        Calculations for your submitted molecule, """+str(identifier)+""", have finished. You may now view the computed properties for this molecule at https://d3tales.as.uky.edu/database/"""+str(identifier)+"""
        
        Thank you, 
        D3TaLES Computational Team"""

        mail = Mail()
        mail.send([email_address], "D3TaLES Molecule Calculations", content)
        return FWAction(update_spec={"identifier": identifier})
