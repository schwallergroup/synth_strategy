from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


HUISGEN_REACTIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "{Huisgen_Cu-catalyzed_1,4-subst}",
    "{Huisgen_Ru-catalyzed_1,5_subst}",
    "{Huisgen_disubst-alkyne}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route employs a Huisgen 1,3-dipolar cycloaddition reaction. This check is based on a predefined list of named reactions including variations for alkynes and alkenes, and different catalysts (e.g., Cu, Ru).
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    click_chemistry_detected = False

    # max_depth would be calculated here based on the route structure.
    # For this example, we pass a placeholder.
    max_depth = 0

    def dfs_traverse(node, depth, max_depth):
        nonlocal click_chemistry_detected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a Huisgen cycloaddition reaction
                is_huisgen_reaction_found = False
                for reaction_name in HUISGEN_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_huisgen_reaction_found = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break # Found one, no need to check others for this rsmi

                if is_huisgen_reaction_found:
                    print(f"Click chemistry detected: Huisgen cycloaddition reaction")
                    click_chemistry_detected = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth, max_depth)

    dfs_traverse(route, 1, max_depth)
    return click_chemistry_detected, findings_json
