import streamlit as st

from mol2png.resolver import fetch_smiles, ResolveError
from mol2png.draw import mol_to_png_bytes
from mol2png.info import fetch_compound_info

import py3Dmol
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem

st.set_page_config(page_title="Mol2PNG", page_icon="🧪", layout="wide")
st.title("🔬 Mol2PNG – Vista 2D + 3D de moléculas + datos")

identifier = st.text_input("Identificador químico (CAS, ChEMBL o SMILES)")
bg_choice  = st.radio("Fondo de la imagen", ("Blanco", "Transparente"), horizontal=True)

if st.button("Generar") and identifier:
    try:
        # 1️⃣ Resuelve SMILES y genera PNG
        smiles = fetch_smiles(identifier.strip())
        img    = mol_to_png_bytes(smiles, transparent=(bg_choice == "Transparente"))
        # 2️⃣ Obtiene datos del compuesto
        info   = fetch_compound_info(smiles)
        # 3️⃣ Genera coordenadas 3D y SDF
        mol3d = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol3d, AllChem.ETKDG())
        sdf_block = Chem.MolToMolBlock(mol3d)
        # 4️⃣ Prepara viewer 3D
        view = py3Dmol.view(width=400, height=400)
        view.addModel(sdf_block, "sdf")
        view.setStyle({"stick": {}})
        view.zoomTo()
        html3d = view._make_html()

        # Layout en dos columnas
        col_img, col_tbl = st.columns((1, 1))

        with col_img:
            st.subheader("🔎 Vista 2D")
            st.image(img, caption=f"{identifier} → {smiles}", use_container_width=False)
            st.download_button("📥 Descargar PNG", img, "mol.png", "image/png")

            st.subheader("🧭 Vista 3D interactiva")
            components.html(html3d, height=420)

        with col_tbl:
            st.subheader("📊 Datos del compuesto")
            if info:
                filas = [{"Propiedad": k, "Valor": v} for k, v in info.items()]
                st.table(filas)
            else:
                st.warning("No se encontraron datos adicionales en PubChem.")

    except ResolveError as e:
        st.error(f"⛔ {e}")
    except Exception as e:  # noqa: BLE001
        st.error(f"⚠️ Error inesperado: {e}")
