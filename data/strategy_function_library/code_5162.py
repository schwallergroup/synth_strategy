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


HETEROCYCLES_OF_INTEREST = ["thiazole", "pyrazole", "pyridine"]

def main(route) -> Tuple[bool, Dict]:
    """Detects a late-stage convergent coupling reaction that joins fragments containing at least two different specified heterocycles from a predefined list. The heterocycles of interest are: thiazole, pyrazole, and pyridine."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, findings_json
        if is_convergent:
            return

        if node["type"] == "reaction" and depth <= 1:
            # Structural Constraint: positional (late_stage)
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "any_reaction",
                    "position": "late_stage",
                    "max_depth": 1
                }
            })

            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            if len(reactants) >= 2:
                # Structural Constraint: count (reactants_in_reaction)
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "reactants_in_reaction",
                        "operator": ">=",
                        "value": 2
                    }
                })

                heterocycles_in_reactants = set()
                for reactant in reactants:
                    for h in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(h, reactant):
                            heterocycles_in_reactants.add(h)
                            # Atomic Check: ring_systems
                            if h not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(h)

                if len(heterocycles_in_reactants) >= 2:
                    # Structural Constraint: count (unique_heterocycle_types_in_reactants)
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "unique_heterocycle_types_in_reactants",
                            "operator": ">=",
                            "value": 2,
                            "from_list": [
                                "thiazole",
                                "pyrazole",
                                "pyridine"
                            ]
                        }
                    })
                    is_convergent = True
                    return

        for child in node.get("children", []):
            if is_convergent:
                break
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # This means it's a chemical node
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return is_convergent, findings_json