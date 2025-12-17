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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving thioether formation
    via mesylate displacement by a thiol.
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

    has_thioether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_thioether_formation, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for S-alkylation reactions
                if checker.check_reaction("S-alkylation of thiols", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("S-alkylation of thiols")
                    print(f"Found S-alkylation reaction: {rsmi}")

                    # Check for mesylate in reactants
                    has_mesylate = False
                    for r in reactants:
                        if checker.check_fg("Mesylate", r):
                            has_mesylate = True
                            findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
                            break
                    print(f"Has mesylate in reactants: {has_mesylate}")

                    # Check for thiol in reactants
                    has_thiol = False
                    for r in reactants:
                        if checker.check_fg("Aliphatic thiol", r):
                            has_thiol = True
                            findings_json["atomic_checks"]["functional_groups"].append("Aliphatic thiol")
                            break
                        if checker.check_fg("Aromatic thiol", r):
                            has_thiol = True
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic thiol")
                            break
                    print(f"Has thiol in reactants: {has_thiol}")

                    # Check for thioether (monosulfide) in product
                    has_thioether = False
                    if checker.check_fg("Monosulfide", product):
                        has_thioether = True
                        findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")
                    print(f"Has thioether in product: {has_thioether}")

                    # If all conditions are met, we've found our strategy
                    if has_mesylate and has_thiol and has_thioether:
                        has_thioether_formation = True
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "S-alkylation of thiols",
                                    "Mesylate",
                                    "Aliphatic thiol",
                                    "Aromatic thiol",
                                    "Monosulfide"
                                ],
                                "description": "A single reaction step must be an 'S-alkylation of thiols' where the reactants contain a Mesylate group and a Thiol group (aliphatic or aromatic), and the product contains a Monosulfide group."
                            }
                        })
                        print(f"Found thioether formation via mesylate displacement: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Thioether formation via mesylate strategy detected: {has_thioether_formation}")
    return has_thioether_formation, findings_json
