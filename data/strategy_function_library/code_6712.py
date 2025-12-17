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

HETEROCYCLES_OF_INTEREST = ["indole", "isoxazole"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy that preserves specific heterocyclic scaffolds, as defined in HETEROCYCLES_OF_INTEREST,
    throughout the entire synthesis. The strategy is flagged if a listed heterocycle is present in the final product
    and is not formed in any preceding synthetic step.
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

    # Initialize tracking variables
    heterocycles_to_track = HETEROCYCLES_OF_INTEREST
    heterocycle_presence = {cycle: False for cycle in heterocycles_to_track}
    preserves_heterocycles = True

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_presence, preserves_heterocycles, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for heterocycles in molecule
            for cycle in heterocycles_to_track:
                if checker.check_ring(cycle, mol_smiles):
                    if cycle not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(cycle)
                    # For the final product (depth 0), record which heterocycles are present
                    if depth == 0:
                        heterocycle_presence[cycle] = True

        elif node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check which heterocycles are present in the product
                product_heterocycles = {}
                for cycle in heterocycles_to_track:
                    if heterocycle_presence[cycle]:  # Only check cycles found in final product
                        product_heterocycles[cycle] = checker.check_ring(cycle, product_smiles)

                # Check if any tracked heterocycle is present in product but missing in all reactants
                for cycle, present_in_product in product_heterocycles.items():
                    if present_in_product:
                        # The heterocycle is in the product, check if it's in any reactant
                        present_in_reactants = False
                        for reactant in reactants_smiles:
                            if checker.check_ring(cycle, reactant):
                                present_in_reactants = True
                                break

                        # If heterocycle is in product but not in any reactant, it's not preserved
                        if not present_in_reactants:
                            preserves_heterocycles = False
                            # If a heterocycle is formed, record the 'ring_formation' atomic check
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception:
                pass

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least one heterocycle is present and preserved
    any_heterocycle_present = any(heterocycle_presence.values())
    strategy_detected = any_heterocycle_present and preserves_heterocycles

    if any_heterocycle_present:
        # Add positional constraint if any heterocycle is present in the final product
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": list(set(c for c, p in heterocycle_presence.items() if p)), # Only add those actually present
                "position": "last_stage",
                "quantifier": "any",
                "condition": "presence"
            }
        })

    if not preserves_heterocycles:
        # Add negation constraint if any heterocycle was formed (not preserved)
        # This implies 'ring_formation' was detected for one of the target heterocycles
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_formation",
                "scope": list(set(c for c, p in heterocycle_presence.items() if p)) # Scope to those that were relevant
            }
        })

    return strategy_detected, findings_json
