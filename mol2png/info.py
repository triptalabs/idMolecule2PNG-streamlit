import requests
from mol2png.config import HTTP_TIMEOUT

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPS = "MolecularFormula,MolecularWeight,IUPACName,IsomericSMILES,InChIKey"

def fetch_compound_info(smiles: str) -> dict[str, str]:
    """Devuelve info de PubChem: fórmula, peso, IUPAC, SMILES e InChIKey."""
    url = f"{PUBCHEM}/compound/smiles/{smiles}/property/{PROPS}/JSON"
    r = requests.get(url, timeout=HTTP_TIMEOUT)
    if r.status_code != 200:
        return {}
    p = r.json()["PropertyTable"]["Properties"][0]

    # Convertir peso a float, si falla usar '-'
    raw_w = p.get("MolecularWeight", None)
    try:
        peso = f"{float(raw_w):.2f}"
    except Exception:
        peso = "-"

    return {
        "Fórmula": p.get("MolecularFormula", "-"),
        "Peso (g/mol)": peso,
        "IUPAC": p.get("IUPACName", "-"),
        "SMILES": p.get("IsomericSMILES", "-"),
        "InChIKey": p.get("InChIKey", "-"),
    }
