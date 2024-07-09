import uuid
import json
import pandas as pd


class ParseExcelMixin:
    """
    Class to process raw UV-Vis data in an Excel format
    Copyright 2021, University of Kentucky
    """
    filepath: str
    data_df: pd.DataFrame
    string_data: pd.DataFrame

    def parse(self):
        """
        Use Pandas to parse the raw data file
        """
        try:
            df = pd.read_excel(self.filepath, header=None, names=['col1', 'col2'])
        except:
            df = pd.read_csv(self.filepath, header=None, names=['col1', 'col2'])

        self.data_df = df.iloc[4:, :].astype(float, errors='ignore').rename(columns={'col1': 'wavelength',
                                                                                     'col2': 'absorbance'})
        self.string_data = df.iloc[:3, :]

    @property
    def integration_time(self):
        """Integration time"""
        query = self.string_data[self.string_data["col1"].str.contains('Integration Time')]['col2'].values
        return query[0] if query else None

    @property
    def date_recorded(self):
        """Data Recorded"""
        query = self.string_data[self.string_data["col1"].str.contains('Timestamp')]['col2'].values
        return query[0] if query else ''

    @property
    def absorbance_data(self):
        """Absorbance data (dict)"""
        return self.data_df.to_dict('list')


class ProcessUvVis(ParseExcelMixin):
    """
    Class to process UV-Vis data files.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, mol_id, metadata=None):
        """
        :param filepath: str, filepath to data file
        :param mol_id: str, molecule ID
        :param metadata: dict, dictionary containing any metadata for this molecule, e.g., {"solvent": "acetonitrile"}
        """
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())
        self.filepath = filepath

        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.solvent = metadata.get("solvent", '')

        self.parse()

    @property
    def no_sql_data(self):
        """
        UV-Vis information in a dictionary that matches the No-SQL schema
        """
        all_data_dict = {
            "date_recorded": self.date_recorded,
            "solvent": self.solvent,
            "instrument": self.instrument,
            "integration_time": self.integration_time,
            "absorbance_data": self.absorbance_data,
        }
        json_data = json.dumps(all_data_dict, default=str)
        return json.loads(json_data)

    @property
    def sql_absorbance_data(self):
        """
        UV-Vis information in a dictionary that matches the SQL AbsorbanceData Table schema
        """
        data = self.absorbance_data
        return [{"uvvis_id": self.uuid,
                 "mol_id": self.mol_id,
                 "wavelength": wavelength,
                 "absorbance": absorbance}
                for wavelength, absorbance in zip(data["wavelength"], data["absorbance"])]

    @property
    def sql_data(self):
        """
        UV-Vis information in a dictionary that matches the SQL UVVisData Table schema
        """
        data_dict = {
            "uvvis_id": self.uuid,
            "mol_id": self.mol_id,
            "date_recorded": self.date_recorded,
            "solvent": self.solvent,
            "instrument": self.instrument,
            "integration_time": self.integration_time,
            "absorbance_data": self.absorbance_data,
        }
        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)
