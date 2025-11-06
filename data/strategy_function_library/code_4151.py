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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route preserves lactam structures throughout the synthesis.
    This means lactam rings are present in the final product and are not modified during synthesis.
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

    # Check if the final product contains a lactam
    final_product_smiles = route.get("smiles", "")
    final_product_has_lactam = False

    if final_product_smiles:
        # Check if the final product contains a lactam using checker functions
        if checker.check_ring("pyrrolidone", final_product_smiles):
            final_product_has_lactam = True
            findings_json["atomic_checks"]["ring_systems"].append("pyrrolidone")
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "pyrrolidone", "position": "last_stage"}})
            print(f"Final product contains lactam: {final_product_smiles}")

    # Track if lactam is preserved throughout synthesis
    lactam_preserved = final_product_has_lactam

    def dfs_traverse(node, depth=0):
        nonlocal lactam_preserved, findings_json

        # For reaction nodes, check if lactam is preserved
        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Extract product
            try:
                product_str = rsmi.split(">")[-1]

                if product_str:
                    # Check if the product contains a lactam using checker functions
                    product_has_lactam = checker.check_ring(
                        "pyrrolidone", product_str
                    )

                    # If the product should have a lactam but doesn't, lactam is not preserved
                    if final_product_has_lactam and not product_has_lactam:
                        lactam_preserved = False
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "absence_of_pyrrolidone_in_intermediate_product"}})
                        print(f"Lactam not preserved at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction for lactam preservation: {e}")

        # Recursively traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return lactam_preserved, findings_json
