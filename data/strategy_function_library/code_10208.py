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


# HETEROCYCLES_OF_INTEREST is a new module-level constant created via Refactoring for Enumeration
HETEROCYCLES_OF_INTEREST = ['benzofuran', 'pyridine']

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a specific set of heterocycles are all present
    in the final product, and at least one of them was formed during the synthesis.
    The heterocycles checked are defined in HETEROCYCLES_OF_INTEREST: ['benzofuran', 'pyridine'].
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

    # State tracking using dictionaries keyed by the heterocycle name
    has_in_product = {h: False for h in HETEROCYCLES_OF_INTEREST}
    is_formed = {h: False for h in HETEROCYCLES_OF_INTEREST}
    formed_heterocycles_count = 0

    def check_heterocycle_formation(reaction_node):
        """Check if a reaction forms any of the heterocycles of interest."""
        nonlocal is_formed, formed_heterocycles_count, findings_json
        try:
            rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            for h in HETEROCYCLES_OF_INTEREST:
                # Check if heterocycle is formed in this reaction
                if checker.check_ring(h, product_smiles) and not any(
                    checker.check_ring(h, r) for r in reactants_smiles
                ):
                    if not is_formed[h]: # Only increment and add to findings if not already marked
                        is_formed[h] = True
                        formed_heterocycles_count += 1
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        except (KeyError, IndexError):
            # Handle cases where rsmi might be missing or malformed
            pass

    def dfs_traverse(node, depth=0):
        """Traverse the synthesis route tree to check for heterocycles."""
        nonlocal has_in_product, findings_json

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if not smiles:
                return

            # Check if this is the final product (root of the tree)
            if depth == 0:
                for h in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(h, smiles):
                        has_in_product[h] = True
                        if h not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(h)

        elif node["type"] == "reaction":
            check_heterocycle_formation(node)

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Depth increases only when NOT traversing from a reaction node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root of the synthesis tree
    dfs_traverse(route)

    # The strategy is detected if:
    # 1. All specified heterocycles are in the final product.
    # 2. At least one of them was formed during the synthesis.
    all_present_in_product = all(has_in_product.values())
    any_formed_in_synthesis = any(is_formed.values())

    result = all_present_in_product and any_formed_in_synthesis

    # Add structural constraints to findings_json if conditions are met
    if has_in_product["benzofuran"]:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "benzofuran", "position": "last_stage"}})
    if has_in_product["pyridine"]:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "pyridine", "position": "last_stage"}})
    
    if formed_heterocycles_count >= 1:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "ring_formation", "operator": ">=", "value": 1}})

    return result, findings_json
