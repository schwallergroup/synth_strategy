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
    This function detects purine scaffold construction in early stage synthesis.
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

    early_stage_construction = False

    def dfs_traverse(node, depth=0):
        nonlocal early_stage_construction, findings_json

        if node["type"] == "reaction" and depth >= 2:  # Early stage (depth >= 2)
            try:
                # Get reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Extract product and reactants
                product_smiles = rsmi.split(">")[-1]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check if product contains purine scaffold
                has_purine_in_product = checker.check_ring("purine", product_smiles)
                if has_purine_in_product:
                    if "purine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("purine")
                print(f"Product contains purine: {has_purine_in_product}")

                # Check if any reactant contains purine scaffold
                reactants_have_purine = False
                for reactant in reactants_smiles:
                    if checker.check_ring("purine", reactant):
                        reactants_have_purine = True
                        if "purine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("purine")
                        print(f"Reactant contains purine: {reactant}")
                        break

                # If product has purine but reactants don't, it's a purine formation
                if has_purine_in_product and not reactants_have_purine:
                    print(f"Detected purine scaffold construction at depth {depth}")
                    early_stage_construction = True
                    # Add structural constraint finding
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "purine_ring_formation",
                            "position": "early_stage"
                        }
                    })
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Early stage purine construction detected: {early_stage_construction}")
    return early_stage_construction, findings_json
