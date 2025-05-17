import re
import requests
from .config import HTTP_TIMEOUT

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL  = "https://www.ebi.ac.uk/chembl/api/data/molecule/{}?format=json"

CAS_RX    = re.compile(r"^\d{2,7}-\d{2}-\d$")
SMILES_RX = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+$")


class ResolveError(Exception):
    """Error al traducir identificadores a SMILES."""


def fetch_smiles(identifier: str) -> str:
    """Devuelve SMILES canónico a partir de CAS, ChEMBL o SMILES."""
    # CAS → PubChem
    if CAS_RX.fullmatch(identifier):
        url = f"{PUBCHEM}/compound/name/{identifier}/property/IsomericSMILES/JSON"
        r = requests.get(url, timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]

    # ChEMBL
    if identifier.lower().startswith("chembl"):
        r = requests.get(CHEMBL.format(identifier.upper()), timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json()["molecule_structures"]["canonical_smiles"]

    # SMILES directo
    if SMILES_RX.fullmatch(identifier):
        return identifier

    raise ResolveError("Identificador no soportado o mal formado.")
