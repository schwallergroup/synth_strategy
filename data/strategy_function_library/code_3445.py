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
    """This function detects a synthetic strategy where a potential benzyl carbamate group is maintained throughout the synthesis. It approximates the presence of a benzyl carbamate by checking for the co-occurrence of a 'Carbamic ester' functional group and a 'benzene' ring in the same molecule at each step."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track presence of benzyl carbamate in each step
    protection_steps = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal protection_steps, max_depth, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check product for benzyl carbamate
            if product_smiles:
                has_cbz_fg = checker.check_fg("Carbamic ester", product_smiles)
                has_cbz_ring = checker.check_ring("benzene", product_smiles)
                has_cbz = has_cbz_fg and has_cbz_ring

                if has_cbz_fg and "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                if has_cbz_ring and "benzene" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("benzene")

                if has_cbz:
                    protection_steps[depth] = True
                    print(f"Found benzyl carbamate at depth {depth} in product: {product_smiles}")

                else:
                    protection_steps[depth] = False

            # Check reactants for benzyl carbamate
            for reactant in reactants:
                if reactant:
                    has_cbz_fg_reactant = checker.check_fg("Carbamic ester", reactant)
                    has_cbz_ring_reactant = checker.check_ring("benzene", reactant)
                    
                    if has_cbz_fg_reactant and "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                    if has_cbz_ring_reactant and "benzene" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("benzene")

                    if has_cbz_fg_reactant and has_cbz_ring_reactant:
                        # If not already marked as protected at this depth
                        if depth not in protection_steps:
                            protection_steps[depth] = True
                            print(f"Found benzyl carbamate at depth {depth} in reactant: {reactant}")

        elif node["type"] == "mol" and not node.get("in_stock", False):
            # Check intermediate molecules
            mol_smiles = node.get("smiles", "")
            if mol_smiles:
                has_cbz_fg_mol = checker.check_fg("Carbamic ester", mol_smiles)
                has_cbz_ring_mol = checker.check_ring("benzene", mol_smiles)

                if has_cbz_fg_mol and "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                if has_cbz_ring_mol and "benzene" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("benzene")

                if has_cbz_fg_mol and has_cbz_ring_mol:
                    protection_steps[depth] = True
                    print(f"Found benzyl carbamate at depth {depth} in molecule: {mol_smiles}")

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if benzyl carbamate protection is maintained throughout
    # We need protection in at least 2 different depths
    print(f"Protection steps: {protection_steps}")
    print(f"Max depth: {max_depth}")

    result = False

    # Strategy criteria:
    # 1. Protection present in at least 2 steps
    # 2. Protection maintained throughout (no gaps)
    if len(protection_steps) < 2:
        print("Less than 2 protection steps found")
        return result, findings_json

    # Check if protection is maintained without gaps
    # Find the first and last depth with protection
    protected_depths = sorted([d for d, v in protection_steps.items() if v])
    if not protected_depths:
        print("No protected depths found")
        return result, findings_json

    first_protected = min(protected_depths)
    last_protected = max(protected_depths)

    # Check if there are any gaps in protection
    for d in range(first_protected, last_protected + 1):
        if d not in protection_steps or not protection_steps[d]:
            print(f"Protection gap found at depth {d}")
            return result, findings_json

    # Ensure protection spans multiple steps (at least 2 depths)
    if last_protected - first_protected < 1:
        print("Protection doesn't span multiple steps")
        return result, findings_json

    print(
        f"Benzyl carbamate protection maintained from depth {first_protected} to {last_protected}"
    )
    result = True

    # Add structural constraints to findings_json if all conditions are met
    if result:
        # Co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Carbamic ester",
                    "benzene"
                ],
                "scope": "same_molecule"
            }
        })
        # Count constraint
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "co-occurrence(Carbamic ester, benzene)",
                "operator": ">=",
                "value": 2,
                "scope": "steps"
            }
        })
        # Sequence constraint
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "event": "co-occurrence(Carbamic ester, benzene)",
                "constraint": "contiguous_block"
            }
        })

    return result, findings_json
