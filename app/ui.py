import streamlit as st
from .resolver import fetch_smiles, ResolveError
from .draw import mol_to_png_bytes

st.set_page_config(page_title="Mol2PNG", page_icon="🧪", layout="centered")
st.title("🔬 Mol2PNG – Vista rápida de moléculas")

identifier = st.text_input("Identificador (CAS, ChEMBL o SMILES)")
bg_choice  = st.radio("Fondo", ("Blanco", "Transparente"), horizontal=True)

if st.button("Generar PNG") and identifier:
    try:
        smiles = fetch_smiles(identifier.strip())
        img    = mol_to_png_bytes(smiles, transparent=(bg_choice == "Transparente"))
        st.image(img, caption=f"{identifier} → {smiles}", use_column_width=False)
        st.download_button("Descargar mol.png", img, "mol.png", "image/png")
    except ResolveError as e:
        st.error(f"⛔ {e}")
    except Exception as e:  # noqa: BLE001
        st.error(f"⚠️ No se pudo generar la imagen: {e}")
