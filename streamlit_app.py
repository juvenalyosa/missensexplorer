# Entry point for `streamlit run streamlit_app.py` and Streamlit Community Cloud.
import runpy
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent / "src"))
runpy.run_path(str(Path(__file__).parent / "src" / "missensexplorer" / "webapp.py"), run_name="__main__")
