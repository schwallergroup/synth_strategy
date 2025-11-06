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
    This function detects a synthetic strategy involving late-stage introduction
    of a trifluoromethyl-containing sulfonamide group.
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

    has_late_stage_cf3_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_cf3_sulfonamide, findings_json

        # Check reaction nodes at late stage (final or penultimate step)
        if node["type"] == "reaction" and depth <= 1:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    print(f"Examining reaction at depth {depth}: {rsmi}")

                    # Properly split reaction SMILES
                    parts = rsmi.split(">")
                    reactants_str = parts[0]
                    product_str = parts[-1]

                    # Skip if product is empty
                    if not product_str:
                        print("Empty product string, skipping")
                        return

                    reactants = reactants_str.split(".")
                    product_mol = Chem.MolFromSmiles(product_str)

                    # Check if this is a sulfonamide formation reaction
                    is_sulfonamide_reaction = False
                    if checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    ):
                        is_sulfonamide_reaction = True
                        if "Sulfonamide synthesis (Schotten-Baumann) primary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) primary amine")
                    
                    if checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    ):
                        is_sulfonamide_reaction = True
                        if "Sulfonamide synthesis (Schotten-Baumann) secondary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) secondary amine")

                    # Check for sulfonyl chloride in reactants and S(=O)(=O)N pattern in product
                    sulfonyl_halide_in_reactants = False
                    for r in reactants:
                        if r and checker.check_fg("Sulfonyl halide", r):
                            sulfonyl_halide_in_reactants = True
                            if "Sulfonyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Sulfonyl halide")
                            break

                    sulfonamide_in_product = False
                    if product_mol is not None and product_mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)N")):
                        sulfonamide_in_product = True
                        if "sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("sulfonamide")

                    if sulfonyl_halide_in_reactants and sulfonamide_in_product:
                        is_sulfonamide_reaction = True

                    print(f"Is sulfonamide reaction: {is_sulfonamide_reaction}")

                    # Check for trifluoromethyl group in reactants
                    has_cf3_reactant = False
                    for r in reactants:
                        if r and checker.check_fg("Trifluoro group", r):
                            has_cf3_reactant = True
                            if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                            break
                    print(f"Has CF3 in reactants: {has_cf3_reactant}")

                    # Check if product contains trifluoromethyl group
                    has_cf3_product = checker.check_fg("Trifluoro group", product_str)
                    if has_cf3_product:
                        if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                    print(f"Has CF3 in product: {has_cf3_product}")

                    # Determine if this is a late-stage trifluoromethyl sulfonamide formation
                    if (
                        is_sulfonamide_reaction
                        and has_cf3_reactant
                        and has_cf3_product
                    ):
                        print(
                            f"Found late-stage trifluoromethyl sulfonamide formation at depth {depth}"
                        )
                        has_late_stage_cf3_sulfonamide = True
                        # Add the structural constraint if all conditions are met
                        if {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "sulfonamide_formation",
                                    "reactant_with_trifluoro_group",
                                    "late_stage_position"
                                ],
                                "description": "A single reaction step must be a sulfonamide formation, have a reactant containing a trifluoro group, and occur at a late stage (depth <= 1)."
                            }
                        } not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "sulfonamide_formation",
                                        "reactant_with_trifluoro_group",
                                        "late_stage_position"
                                    ],
                                    "description": "A single reaction step must be a sulfonamide formation, have a reactant containing a trifluoro group, and occur at a late stage (depth <= 1)."
                                }
                            })
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_cf3_sulfonamide, findings_json
