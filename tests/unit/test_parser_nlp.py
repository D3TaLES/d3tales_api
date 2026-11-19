import json
import sys
import types
from datetime import datetime
from types import SimpleNamespace

import pandas as pd
import pytest

# NLP integrations are optional in the runtime environment; these tests exercise
# parsing and transformation logic without downloading articles or using Elsevier.
_chem = types.ModuleType("chemdataextractor2")
_chem.Document = object
_chem_text = types.ModuleType("chemdataextractor2.doc.text")
_chem_text.Paragraph = type("Paragraph", (), {})
_chem_doc = types.ModuleType("chemdataextractor2.doc")
_chem_doc.text = _chem_text
_elsdoc = types.ModuleType("elsapy.elsdoc")
_elsdoc.AbsDoc = object
_elsdoc.Paragraph = type("Paragraph", (), {})
_elsclient = types.ModuleType("elsapy.elsclient")
_elsclient.ElsClient = object
_els_package = types.ModuleType("elsapy")
_els_package.elsdoc = _elsdoc
_els_package.elsclient = _elsclient
_article_module = types.ModuleType("articledownloader.articledownloader")
_article_module.ArticleDownloader = object
_article_package = types.ModuleType("articledownloader")
_article_package.articledownloader = _article_module
sys.modules.setdefault("chemdataextractor2", _chem)
sys.modules.setdefault("chemdataextractor2.doc", _chem_doc)
sys.modules.setdefault("chemdataextractor2.doc.text", _chem_text)
sys.modules.setdefault("elsapy", _els_package)
sys.modules.setdefault("elsapy.elsdoc", _elsdoc)
sys.modules.setdefault("elsapy.elsclient", _elsclient)
sys.modules.setdefault("articledownloader", _article_package)
sys.modules.setdefault("articledownloader.articledownloader", _article_module)

from d3tales_api.Processors import parser_nlp
from d3tales_api.Processors.parser_nlp import ProcessNlp


class CapturedNlp:
    def __init__(self, instance=None, **kwargs):
        self.instance = instance
        self.kwargs = kwargs


@pytest.mark.unit
def test_process_nlp_requires_doi_and_builds_data_dict():
    with pytest.raises(ValueError, match="requires a DOI"):
        ProcessNlp.__new__(ProcessNlp).__init__()

    obj = ProcessNlp.__new__(ProcessNlp)
    obj.doi = "10/example"
    obj.basic_info = {"publisher": "Publisher"}
    obj.instance = {"title": "Article"}
    obj.article_download = True
    obj.main_text = {"text": ["paragraph"]}
    obj.download_path = "article.html"
    assert obj.data_dict == {
        "publisher": "Publisher",
        "title": "Article",
        "_id": "10/example",
        "main_text": {"text": ["paragraph"]},
        "pdf_location": obj.data_dict["pdf_location"],
    }


@pytest.mark.unit
def test_author_affiliation_shapes():
    assert ProcessNlp.get_affiliations([{"@id": "1"}, {"@id": "2"}]) == ["1", "2"]
    assert ProcessNlp.get_affiliations({"@id": "3"}) == ["3"]
    assert ProcessNlp.get_affiliations("unknown") == []


@pytest.mark.unit
def test_basic_info_uses_elsevier_record_and_maps_authors(monkeypatch):
    class Schema:
        def Nlp(self):
            result = SimpleNamespace()
            result.as_dict = lambda: {key: value for key, value in vars(result).items() if key != "as_dict"}
            return result

        def Author(self):
            return SimpleNamespace(from_json=lambda value: json.loads(value))

    article = SimpleNamespace(
        title="Title",
        data={
            "coredata": {
                "dc:publisher": "Publisher",
                "prism:publicationName": "Journal",
                "prism:coverDate": "2025",
                "prism:doi": "10/example",
                "dc:description": "Abstract",
            },
            "authors": {"author": [{"@auid": "a1", "ce:surname": "Smith", "affiliation": {"@id": "org"}}]},
        },
        read=lambda client: True,
    )
    monkeypatch.setattr(parser_nlp, "Schema2Class", lambda **kwargs: Schema())
    monkeypatch.setattr(parser_nlp, "ElsClient", lambda key: ("client", key))
    monkeypatch.setattr(parser_nlp, "AbsDoc", lambda uri: article)
    obj = ProcessNlp.__new__(ProcessNlp)
    obj.els_apikey = "key"
    obj.doi = "10/example"
    obj.s2c = Schema()

    result = obj.get_basic_info()

    assert result["publisher"] == "Publisher"
    assert result["journal"] == "Journal"
    assert result["title"] == "Title"
    assert result["authors"][0]["affiliations"] == ["org"]
    assert result["_id"] == "10/example"


@pytest.mark.unit
def test_basic_info_raises_when_abstract_is_unavailable(monkeypatch):
    class Schema:
        def Nlp(self):
            return SimpleNamespace()

        def Author(self):
            return SimpleNamespace()

    article = SimpleNamespace(read=lambda client: False)
    monkeypatch.setattr(parser_nlp, "Schema2Class", lambda **kwargs: Schema())
    monkeypatch.setattr(parser_nlp, "ElsClient", lambda key: None)
    monkeypatch.setattr(parser_nlp, "AbsDoc", lambda uri: article)
    obj = ProcessNlp.__new__(ProcessNlp)
    obj.els_apikey = None
    obj.doi = "10/example"
    obj.s2c = Schema()
    with pytest.raises(ConnectionError, match="Not able to retrieve"):
        obj.get_basic_info()


@pytest.mark.unit
def test_main_text_extracts_paragraphs_before_references(monkeypatch):
    class FakeParagraph:
        def __init__(self, text):
            self.text = text

    elements = [
        SimpleNamespace(text="Heading"),
        SimpleNamespace(text="Abstract"),
        FakeParagraph("First paragraph"),
        FakeParagraph("Second paragraph"),
        SimpleNamespace(text="References"),
        FakeParagraph("Not included"),
    ]

    class Data:
        def __init__(self):
            self.text = []

        def as_dict(self):
            return {"text": self.text}

    monkeypatch.setattr(parser_nlp, "Document", SimpleNamespace(from_file=lambda path: SimpleNamespace(elements=elements)))
    monkeypatch.setattr(parser_nlp, "Paragraph", FakeParagraph)
    obj = ProcessNlp.__new__(ProcessNlp)
    obj.download_path = "paper.html"
    obj.doi = "10/example"
    obj.s2c = SimpleNamespace(Data=Data)
    assert obj.get_main_text() == {"text": ["First paragraph", "Second paragraph"]}


@pytest.mark.unit
@pytest.mark.parametrize(
    "publisher, expected_mode",
    [
        ("American Chemical Society", "acs"),
        ("Royal Society of Chemistry", "rsc"),
        ("Elsevier B.V.", "elsevier"),
        ("Nature Research", "nature"),
        ("IOP Publishing Ltd", "ecs"),
        ("Wiley-VCH Verlag", "wiley"),
        ("Korean Electrochemical Society", "ecs"),
        ("Other", "elsevier"),
    ],
)
def test_download_article_selects_publisher_mode(tmp_path, monkeypatch, publisher, expected_mode):
    modes = []

    class Downloader:
        def __init__(self, els_api_key):
            pass

        def get_html_from_doi(self, doi, writefile, mode):
            modes.append(mode)
            writefile.write(b"<html/>")

    monkeypatch.setattr(parser_nlp, "ArticleDownloader", Downloader)
    obj = ProcessNlp.__new__(ProcessNlp)
    obj.els_apikey = "key"
    obj.publisher = publisher
    obj.doi = "10/example"
    obj.download_path = str(tmp_path / "paper.html")
    obj.download_article()
    assert modes == [expected_mode]
    assert (tmp_path / "paper.html").read_bytes() == b"<html/>"


@pytest.mark.unit
def test_dataframe_conversion_formats_properties_and_rejects_bad_columns():
    columns = ["molecule", "property", "value", "unit", "line_number", "parent_sentence", "notes"]
    frame = pd.DataFrame(
        [["ferrocene", "Oxidation Potential", "0.5", "V", 3, "oxidized.", "measured"]],
        columns=columns,
    )
    result = ProcessNlp.from_dataframe.__func__(
        CapturedNlp,
        frame,
        "model",
        doi="https://doi.org/10/example",
        base_instance={"source": "paper"},
        date_generated=datetime(2024, 1, 1),
    )
    record = result.instance
    assert record["_id"] == "10/example"
    prop = record["extracted_molecules"][0]["extracted_properties"]["oxidation_potential"][0]
    assert prop["value"] == pytest.approx(0.5)
    assert prop["conditions"]["nlp_model"] == "model"
    assert prop["line_number"] == 3

    with pytest.raises(SyntaxError, match="unexpected columns"):
        ProcessNlp.from_dataframe.__func__(
            CapturedNlp,
            frame.assign(extra=1),
            "model",
            doi="10/example",
            base_instance={},
        )
    with pytest.raises(ValueError, match="unexpected properties"):
        ProcessNlp.from_dataframe.__func__(
            CapturedNlp,
            frame.assign(property="unknown"),
            "model",
            doi="10/example",
            base_instance={},
        )


@pytest.mark.unit
def test_html_and_property_list_factories_validate_input_and_delegate(monkeypatch):
    captured = []
    monkeypatch.setattr(ProcessNlp, "from_dataframe", classmethod(
        lambda cls, frame, model, **kwargs: captured.append((frame, model, kwargs)) or "created"
    ))
    assert ProcessNlp.from_prop_list(
        [{"molecule": "A", "property": "solubility", "value": "1", "unit": "M",
          "line_number": 1, "parent_sentence": "x", "notes": ""}],
        "model",
        doi="10/example",
    ) == "created"
    assert captured[-1][1] == "model"
    with pytest.raises(ValueError, match="no propreties"):
        ProcessNlp.from_prop_list([], "model", doi="10/example")

    monkeypatch.setattr(parser_nlp.pd, "read_html", lambda html: [pd.DataFrame({"value": [1]})])
    assert ProcessNlp.from_html("<table></table>", "model", doi="10/example") == "created"
    monkeypatch.setattr(parser_nlp.pd, "read_html", lambda html: [])
    with pytest.raises(ValueError, match="no tables"):
        ProcessNlp.from_html("<html/>", "model", doi="10/example")
    monkeypatch.setattr(parser_nlp.pd, "read_html", lambda html: [pd.DataFrame(), pd.DataFrame()])
    with pytest.raises(ValueError, match="2 tables"):
        ProcessNlp.from_html("<html/>", "model", doi="10/example")
