import pathlib
import requests
import pandas as pd
import io
import typing as T
from enum import Enum, auto
import urllib.parse
from itertools import islice
import bioinfdatacollector.query_generator as QueryGenerator

import bioinfdatacollector.rhea.config as RC
from bioinfdatacollector.rhea.query_generator import RheaQueryBuilder
from bioinfdatacollector.sparql_query import SelectQuery
import bioinfdatacollector.uniprot.config as UC
from bioinfdatacollector.uniprot.query_generator import UniprotQueryBuilder


def batched(iterable, n):
    "Batch data into lists of length n. The last batch may be shorter."
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := list(islice(it, n)):
        yield batch


class Repository(Enum):
    RHEA = auto()
    UNIPROT = auto()


def collect_data(query: SelectQuery, url: str) -> pd.DataFrame:
    response = requests.get(
        url, params={"query": query.get_pretty_text(), "format": "csv"}, stream=True
    )
    response.raise_for_status()
    if response.ok:
        data_frame = pd.read_csv(
            io.StringIO(response.text), sep=",", header=0, skipinitialspace=True
        )
    return data_frame


TConfig = T.TypeVar("TConfig")


def get_uniprot(features):
    uniprot_cfg = UC.UniprotSearchConfig(
        data_selector=UC.UniprotSelection([UC.Feature.PROTEIN] + features),
        data_filter=UC.UniprotSearchFilter(pfams=["PF00067"]),
    )

    query = get_query(uniprot_cfg, UniprotQueryBuilder)
    uniprot_data = collect_data(query, UC.URL)
    uniprot_data["protein"] = uniprot_data["protein"].apply(
        lambda p: urllib.parse.urlparse(p).path.split("/")[2]
    )
    uniprot_data = uniprot_data.set_index("protein", drop=True)
    return uniprot_data


def get_query(
    config: TConfig,
    query_builder_ctor: T.Callable[
        [TConfig], QueryGenerator.SparqlQueryBuilder[TConfig]
    ],
) -> SelectQuery:
    builder = query_builder_ctor(config)
    return builder.get_query()


def save_data(data: pd.DataFrame, path: pathlib.Path) -> None:
    with open(path, mode="bw") as writer:
        if path.suffix == ".xlsx":
            data.to_excel(writer, index=True, engine="openpyxl")
        if path.suffix == ".csv":
            data.to_csv(writer, index=True)


if __name__ == "__main__":

    print("Fetching taxa")

    taxa = get_uniprot(
        [
            UC.Feature.ORGANISM,
            UC.Feature.KINGDOM,
            UC.Feature.SUPERKINGDOM,
        ]
    )
    save_data(taxa, pathlib.Path("taxaNAT.csv"))

    print("Fetching names")
    names = get_uniprot([UC.Feature.NAME, UC.Feature.SUBMITTED_NAME])
    save_data(names, pathlib.Path("namesNAT.csv"))
    #cytochrome_names = names[names['full_recommended_name'].str.contains('cytochrome', case=False) | names['full_submitted_name'].str.contains('cytochrome', case=False)]
    #cytochrome_names.to_csv('cytochrome_names.csv', index=False)



    print("Fetching pfam")
    pfams = get_uniprot([UC.Feature.PFAM])
    pfams["pfam"] = pfams["pfam"].apply(
        lambda p: urllib.parse.urlparse(p).path.split("/")[2]
    )
    save_data(pfams, pathlib.Path("pfams.csv"))

    print("Fetching sequences")
    sequences = get_uniprot([UC.Feature.SEQUENCE])
    #I would want to fetch lengths as well but not working
    save_data(sequences, pathlib.Path("sequencesNAT.csv"))

    uniprot_data = get_uniprot([UC.Feature.REACTION])

    uniprot_data["reaction"] = (
        uniprot_data["reaction"]
        .apply(lambda r: urllib.parse.urlparse(r).path[1:])
        .astype(str)
    )
    reaction_ids = uniprot_data["reaction"].unique()

    configs = map(
        lambda rs: RC.RheaSearchConfig(
            RC.RheaSelection(
                [
                    RC.Feature.REACTION,
                    RC.Feature.REACTION_SIDE,
                    RC.Feature.CHEBI,
                    RC.Feature.SMILES,
                ]
            ),
            RC.RheaSearchFilter(reactions=rs),
        ),
        batched(reaction_ids, 20),
    )

    print("Fetching reactions")
    queries = list(map(lambda c: get_query(c, RheaQueryBuilder), configs))
    results = map(lambda q: collect_data(q, RC.URL), queries)

    rhea_data = pd.concat(results)
    rhea_data["reaction"] = (
        rhea_data["reaction"]
        .apply(lambda r: urllib.parse.urlparse(r).path[1:])
        .astype(str)
    )
    rhea_data["cofactor"] = rhea_data["cofactor"].apply(
        lambda c: urllib.parse.urlparse(c).path.split("/")[2]
    )
    rhea_data = rhea_data.set_index("reaction", drop=True)
    rhea_data["reaction_side_order"] = rhea_data["reaction_side_order"].apply(
        lambda s: s[0]
    )
    rhea_data["chebi"] = rhea_data["chebi"].apply(
        lambda c: urllib.parse.urlparse(c).path.split("/")[2]
    )
    data = uniprot_data.join(rhea_data, on="reaction")
    save_data(data, pathlib.Path("reaction_dataNAT.csv"))
