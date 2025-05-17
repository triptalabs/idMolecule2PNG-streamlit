import streamlit as st

from mol2png.resolver import fetch_smiles, ResolveError
from mol2png.draw import mol_to_png_bytes
from mol2png.info import fetch_compound_info

st.set_page_config(page_title="Mol2PNG", page_icon="🧪", layout="wide")
st.title("🔬 Mol2PNG – Vista rápida de moléculas + datos")

identifier = st.text_input("Identificador químico (CAS, ChEMBL o SMILES)")
bg_choice  = st.radio("Fondo de la imagen", ("Blanco", "Transparente"), horizontal=True)

if st.button("Generar") and identifier:
    try:
        smiles = fetch_smiles(identifier.strip())
        img    = mol_to_png_bytes(smiles, transparent=(bg_choice == "Transparente"))
        info   = fetch_compound_info(smiles)

        col_img, col_tbl = st.columns((1, 1))

        with col_img:
            st.image(
                img,
                caption=f"{identifier} → {smiles}",
                use_container_width=False
            )
            st.download_button("📥 Descargar mol.png", img, "mol.png", "image/png")

        with col_tbl:
            st.subheader("Datos del compuesto")
            if info:
                # Reemplazamos el dict por una lista de dicts con headers 'Propiedad' y 'Valor'
                filas = [{"Propiedad": prop, "Valor": val} for prop, val in info.items()]
                st.table(filas)
            else:
                st.warning("No se encontraron datos adicionales en PubChem.")

    except ResolveError as e:
        st.error(f"⛔ {e}")
    except Exception as e:  # noqa: BLE001
        st.error(f"⚠️ Error inesperado: {e}")
