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
    This function detects a synthetic strategy involving sequential functionalization
    of a pyrimidine core.
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

    result = False
    pyrimidine_modifications = []

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json
        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for pyrimidine in product
                if checker.check_ring("pyrimidine", product):
                    findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                    # Check if any reactant contains pyrimidine
                    for reactant in reactants:
                        if checker.check_ring("pyrimidine", reactant):
                            findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                            # If pyrimidine is in both reactant and product, it's a modification
                            print(f"Pyrimidine modification detected at depth {depth}")
                            pyrimidine_modifications.append(depth)
                            findings_json["atomic_checks"]["named_reactions"].append("pyrimidine_modification")
                            break

        # Traverse children with increased depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 2 modifications
    if len(pyrimidine_modifications) < 2:
        print("Less than 2 pyrimidine modifications found")
        result = False
    else:
        # Sort modifications by depth to analyze the sequence
        # In retrosynthetic analysis, we're moving from product to reactants
        # So we're looking at the synthesis in reverse
        pyrimidine_modifications.sort()

        print(f"Found pyrimidine modifications at depths: {pyrimidine_modifications}")

        # Check if we have sequential modifications (not necessarily at adjacent depths)
        # The key is that we have multiple modifications on the pyrimidine core
        # throughout the synthesis route
        result = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "pyrimidine_modification",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
