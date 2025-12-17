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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects pyrimidine ring formation in the middle of the synthesis.
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

    pyrimidine_formation = False
    formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal pyrimidine_formation, formation_depth, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyrimidine
                product_has_pyrimidine = checker.check_ring("pyrimidine", product)

                if product_has_pyrimidine:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                    # Check if reactants don't have pyrimidine
                    reactants_have_pyrimidine = False
                    for reactant in reactants:
                        if checker.check_ring("pyrimidine", reactant):
                            reactants_have_pyrimidine = True
                            break

                    if not reactants_have_pyrimidine:
                        # This is a pyrimidine formation reaction
                        pyrimidine_formation = True
                        formation_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # This means it's a 'chemical' node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if pyrimidine formation occurs in the middle of synthesis (depth 1-4)
    middle_stage = formation_depth is not None and 1 <= formation_depth <= 4

    result = pyrimidine_formation and middle_stage

    if pyrimidine_formation:
        # This implies 'ring_formation:pyrimidine' occurred
        # The 'not_last_stage' constraint is implicitly checked by 'formation_depth is not None'
        # and the '1 <= formation_depth' part of 'middle_stage'
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_formation:pyrimidine",
                "position": "not_last_stage"
            }
        })

    if formation_depth is not None and formation_depth <= 4:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "depth_of_ring_formation:pyrimidine",
                "operator": "<=",
                "value": 4
            }
        })

    return result, findings_json
