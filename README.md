# idMolecule2PNG

> Convierte identificadores químicos (CAS, ChEMBL o SMILES) en vistas 2D y 3D interactivas, con datos detallados del compuesto.
> *Powered by Streamlit, RDKit & py3Dmol.*

![Social preview](./banner.png)

---

## 🚀 Demo en vivo

👉 [https://idmolecule2png-streamlit.streamlit.app](https://idmolecule2png-streamlit.streamlit.app)

---

## ✨ Características

* **Resolución de identificadores** (CAS / ChEMBL / SMILES) → SMILES canónico.
* **Renderizado 2D** con RDKit: PNG 400×400 px (fondo blanco o transparente).
* **Viewer 3D interactivo** con py3Dmol: rotación y zoom aplicados sobre coordenadas generadas por RDKit.
* **Datos del compuesto** en tabla **Propiedad / Valor**: fórmula, peso molecular, nombre IUPAC, InChIKey, SMILES.
* **Descarga inmediata** de la imagen 2D (`mol.png`).
* **Cache de datos** con `st.cache_data` para acelerar consultas repetidas a PubChem.
* **Mínimas dependencias**: Streamlit, rdkit‑pypi, requests, numpy, pillow, py3Dmol.

---

## 📦 Instalación local

```bash
# 1. Clona el repositorio
git clone https://github.com/triptalabs/idMolecule2PNG-streamlit.git
cd idMolecule2PNG-streamlit

# 2. Crea y activa un entorno virtual
python -m venv .venv
# Windows
.venv\Scripts\Activate.ps1
# macOS / Linux
source .venv/bin/activate

# 3. Instala dependencias runtime
pip install --upgrade pip
pip install -r requirements.txt
```

---

## 🏃‍♂️ Ejecución

```bash
streamlit run ui.py
```

Abre [http://localhost:8501](http://localhost:8501) en tu navegador.

---

## 🛠 Desarrollo

Para desarrollo y pruebas:

```bash
pip install -r requirements-dev.txt
pre-commit install
pytest -q
```

---

## 📁 Estructura del proyecto

```
idMolecule2PNG-streamlit/
├─ mol2png/            # Lógica de negocio
│  ├─ __init__.py
│  ├─ config.py        # Parámetros y timeouts
│  ├─ resolver.py      # fetch_smiles()
│  ├─ draw.py          # mol_to_png_bytes()
│  └─ info.py          # fetch_compound_info()
│
├─ ui.py               # Launcher Streamlit (2D + 3D + tabla)
├─ requirements.txt    # Dependencias runtime mínimas
├─ requirements-dev.txt# Dependencias de desarrollo
└─ runtime.txt         # Python 3.11 para Streamlit Cloud
```

---

## 📄 API interna

| Función                                       | Descripción                                         |
| --------------------------------------------- | --------------------------------------------------- |
| `fetch_smiles(identifier: str) -> str`        | Convierte CAS / ChEMBL / SMILES en SMILES canónico. |
| `fetch_compound_info(smiles: str) -> dict`    | Obtiene datos del compuesto (PubChem).              |
| `mol_to_png_bytes(smiles, size, transparent)` | Genera imagen PNG RDKit en memoria.                 |

---

## 🌐 Despliegue

### Streamlit Cloud

1. Push a `main` con `ui.py` en la raíz.
2. Asegúrate de tener `runtime.txt` (`python-3.11`) y `requirements.txt` actualizados.
3. En Streamlit Community Cloud, selecciona **Main file** = `ui.py` y despliega.

### VPS + Nginx

Sigue el **proxy inverso** y **Certbot** descritos en la [Wiki del proyecto](./DEVELOPMENT.md).

---

## 🤝 Contribuciones

Las contribuciones son bienvenidas.

1. Haz *fork* del repo.
2. Crea una rama (`git checkout -b feat/mi-feature`).
3. Commit con convención: `feat: descripción breve`.
4. Push y abre *Pull Request*.

---

## 📝 Licencia

Este proyecto está bajo la licencia [MIT](./LICENSE).
