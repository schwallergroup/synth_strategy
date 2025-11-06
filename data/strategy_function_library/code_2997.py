from typing import Tuple, Dict, List
import copy
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


HALOGEN_FGS = [
    "Aromatic halide",
    "Primary halide",
    "Secondary halide",
    "Tertiary halide",
    "Alkenyl halide",
    "Haloalkyne",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage halogen functionalization strategy. This is defined as a halogen handle being present in an intermediate, carried through at least two synthetic steps, and then consumed in one of the final two steps of the synthesis.
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

    halogen_at_depth = {}
    halogen_used_at_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=1):
        nonlocal halogen_used_at_depth, max_depth, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            for fg in HALOGEN_FGS:
                if checker.check_fg(fg, mol_smiles):
                    halogen_at_depth[depth] = halogen_at_depth.get(depth, False) or True
                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg)

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_has_halogen = False
            for r in reactants:
                for fg in HALOGEN_FGS:
                    if checker.check_fg(fg, r):
                        reactant_has_halogen = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

            product_has_halogen = False
            for fg in HALOGEN_FGS:
                if checker.check_fg(fg, product):
                    product_has_halogen = True
                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg)

            if reactant_has_halogen and not product_has_halogen:
                if halogen_used_at_depth is None or depth < halogen_used_at_depth:
                    halogen_used_at_depth = depth
                    if "halogen_consumption" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("halogen_consumption")

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = False

    if halogen_used_at_depth is None:
        return result, findings_json

    preserved_steps = 0
    for d in range(halogen_used_at_depth + 1, max_depth + 1):
        if halogen_at_depth.get(d, False):
            preserved_steps += 1

    if halogen_used_at_depth <= 2:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "halogen_consumption",
                "position": "last_two_stages"
            }
        })

    if preserved_steps >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "halogen_preservation_step",
                "operator": ">=",
                "value": 2
            }
        })

    result = preserved_steps >= 2 and halogen_used_at_depth <= 2

    return result, findings_json
