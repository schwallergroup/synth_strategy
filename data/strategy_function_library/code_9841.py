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
    Detects if a cyano group is present in the target molecule and also in at least one of its precursor molecules. This acts as a proxy to identify routes where the cyano group is carried over from an earlier stage, rather than being introduced in the final step.
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

    cyano_maintained = False
    target_has_cyano = False

    def dfs_traverse(node, depth=0):
        nonlocal cyano_maintained, target_has_cyano, findings_json

        if node["type"] == "mol":
            if node.get("in_stock", False):
                return

            has_cyano = checker.check_fg("Nitrile", node["smiles"])
            if has_cyano:
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

            if depth == 0 and has_cyano:
                target_has_cyano = True
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "Nitrile",
                        "position": "last_stage"
                    }
                })

            if depth >= 1 and has_cyano:
                cyano_maintained = True
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "Nitrile",
                        "position": "not_last_stage"
                    }
                })

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction":  # If current node is 'mol' (chemical), depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = target_has_cyano and cyano_maintained

    return result, findings_json
