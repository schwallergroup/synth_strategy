"""Minimal Streamlit UI for matching synthesis routes against the strategy function library.

Run via the `synth-strategy-ui` console script, or directly with:
    streamlit run src/synth_strategy/webui.py
"""

import ast
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

from .api import SynthStrategyAPI

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FUNCTIONS_DIR = REPO_ROOT / "data" / "strategy_function_library"


def route_from_reaction_smiles(reaction_smiles_list: List[str]) -> Dict[str, Any]:
    """Build a linear synthesis-route tree from an ordered (target-first) list of reaction SMILES."""
    target = reaction_smiles_list[0].split(">>")[1]
    node = {"type": "mol", "smiles": target, "metadata": {"target_smiles": target}}
    current = node
    for rsmi in reaction_smiles_list:
        if ">>" not in rsmi:
            raise ValueError(f"Reaction SMILES missing '>>' separator: {rsmi!r}")
        reaction_node = {"type": "reaction", "metadata": {"mapped_reaction_smiles": rsmi}, "children": []}
        current["children"] = [reaction_node]
        current = reaction_node
    return node


def extract_main_docstring(function_file: Path) -> str:
    """Extract the docstring of the `main` function from a strategy function source file."""
    try:
        tree = ast.parse(function_file.read_text())
    except (OSError, SyntaxError) as e:
        return f"Error parsing file: {e}"
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "main":
            return ast.get_docstring(node) or "No docstring found."
    return "No docstring found."


def main() -> None:
    import streamlit as st

    st.set_page_config(page_title="SynthStrategy Explorer", layout="wide")
    st.title("SynthStrategy: Strategy Function Matcher")
    st.caption(
        "Paste a linear retrosynthetic route (or a raw route JSON) and check which "
        "strategy functions in the library match it."
    )

    mode = st.radio("Input mode", ["Reaction SMILES (one per line)", "Route JSON"])

    if mode == "Reaction SMILES (one per line)":
        text = st.text_area(
            "Reaction SMILES, ordered from final product back to starting materials",
            placeholder="CC(=O)Cl.NCc1ccccc1>>CC(=O)NCc1ccccc1",
            height=150,
        )
    else:
        text = st.text_area("Route JSON (single route object)", height=300)

    if st.button("Run strategy annotation"):
        if not text or not text.strip():
            st.error("Please provide some input before running.")
            return

        try:
            if mode == "Reaction SMILES (one per line)":
                lines = [line.strip() for line in text.splitlines() if line.strip()]
                if not lines:
                    st.error("No reaction SMILES lines found.")
                    return
                route = route_from_reaction_smiles(lines)
            else:
                route = json.loads(text)
        except (ValueError, json.JSONDecodeError) as e:
            st.error(f"Could not parse input: {e}")
            return

        if not DEFAULT_FUNCTIONS_DIR.is_dir():
            st.error(f"Strategy function library not found at: {DEFAULT_FUNCTIONS_DIR}")
            return

        try:
            with st.spinner("Running strategy annotation..."):
                api = SynthStrategyAPI()
                annotated = api.annotate_strategies(
                    routes=[route], functions_dir=str(DEFAULT_FUNCTIONS_DIR)
                )
        except Exception as e:
            st.error(f"Annotation failed: {e}")
            return

        matched = annotated[0].get("passing_functions", {}) if annotated else {}
        st.subheader(f"Matched {len(matched)} strategy function(s)")
        if not matched:
            st.info("No strategy functions matched this route.")
        for name in matched:
            docstring = extract_main_docstring(DEFAULT_FUNCTIONS_DIR / name)
            with st.expander(name):
                st.markdown(docstring)
                st.json(matched[name])


def launch() -> None:
    """Console-script entrypoint: shells out to `streamlit run` on this module."""
    subprocess.run([sys.executable, "-m", "streamlit", "run", str(Path(__file__).resolve())], check=True)


if __name__ == "__main__":
    main()
