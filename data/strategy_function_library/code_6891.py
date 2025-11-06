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
    This function detects a synthetic strategy involving late-stage
    ester hydrolysis to form a carboxylic acid.
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

    ester_hydrolysis_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_detected, findings_json

        if ester_hydrolysis_detected:
            return

        if (
            node["type"] == "reaction" and depth <= 2
        ): 
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # The checker.check_reaction function is the only reliable method here.
                # The original code's alternative checks were flawed and produced false positives.
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    ester_hydrolysis_detected = True
                    findings_json["atomic_checks"]["named_reactions"].append("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                            "max_depth_from_root": 2
                        }
                    })
                    return

        # Recursively process children with increased depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not a reaction (i.e., chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage ester hydrolysis strategy detected: {ester_hydrolysis_detected}")
    return ester_hydrolysis_detected, findings_json
