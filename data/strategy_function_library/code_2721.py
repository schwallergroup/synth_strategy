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


SULFONAMIDE_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route employs a late-stage sulfonamide formation
    as the final step (depth 0) or penultimate step (depth 1).
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

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            if depth <= 1:  # Final or penultimate step (late-stage)
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains sulfonamide
                product_has_sulfonamide = checker.check_fg("Sulfonamide", product_smiles)
                if product_has_sulfonamide:
                    if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

                if product_has_sulfonamide:
                    # Check if sulfonamide is formed in this step (not present in reactants)
                    sulfonamide_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Sulfonamide", reactant):
                            sulfonamide_in_reactants = True
                            break

                    # Check if this is a recognized sulfonamide formation reaction
                    is_sulfonamide_reaction = False
                    for r in SULFONAMIDE_REACTIONS:
                        if checker.check_reaction(r, rsmi):
                            is_sulfonamide_reaction = True
                            if r not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r)
                            break

                    if not sulfonamide_in_reactants and is_sulfonamide_reaction:
                        result = True
                        # Add structural constraints if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "sulfonamide_formation_step",
                                "position": "last_or_penultimate_stage"
                            }
                        })
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "single_reaction_step",
                                "targets": [
                                    "Sulfonamide synthesis (Schotten-Baumann)",
                                    "functional_group_formation"
                                ]
                            }
                        })

        for child in node.get("children", []):
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return result, findings_json
