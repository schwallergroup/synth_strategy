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
    Detects synthesis routes involving purine scaffold construction.
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

    purine_construction = False

    def dfs_traverse(node, depth=0):
        nonlocal purine_construction, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains purine scaffold
                if checker.check_ring("purine", product_smiles):
                    if "purine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("purine")
                    # Check if reactants don't have the purine scaffold
                    reactants_with_purine = sum(
                        1 for r in reactants_smiles if checker.check_ring("purine", r)
                    )

                    if reactants_with_purine == 0:
                        purine_construction = True
                        # Add structural constraint if purine construction is detected
                        # This corresponds to the 'purine_ring_formation' constraint
                        if {"type": "count", "details": {"target": "purine_ring_formation", "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "purine_ring_formation", "operator": ">=", "value": 1}})

            except Exception as e:
                # Handle any errors in SMILES parsing or structure checking
                pass

        # Determine the new depth based on the current node's type
        new_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return purine_construction, findings_json
