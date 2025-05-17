from io import BytesIO
from typing import Tuple

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def mol_to_png_bytes(smiles: str,
                     size: Tuple[int, int] = (400, 400),
                     transparent: bool = False) -> bytes:
    """Convierte SMILES en PNG bruto."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES inválido")

    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    drawer.drawOptions().clearBackground = not transparent
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
