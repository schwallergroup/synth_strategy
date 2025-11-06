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

HETEROCYCLES_OF_INTEREST = [
    "pyrazole",
    "isoxazole",
    "oxazole",
    "thiazole",
    "imidazole",
    "triazole",
    "tetrazole",
    "furan",
    "pyrrole",
    "thiophene",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the preservation of a specific heterocyclic core throughout a synthesis. This is confirmed if at least two molecules in the synthetic tree, including the final product, contain one of the specified heterocycles: pyrazole, isoxazole, oxazole, thiazole, imidazole, triazole, tetrazole, furan, pyrrole, or thiophene.
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

    # Track heterocycle-containing molecules in the main synthetic pathway
    main_pathway_molecules = []
    heterocycle_molecules = []
    max_depth = 0
    result = False

    # List of heterocycles is now a module-level constant HETEROCYCLES_OF_INTEREST

    def dfs_traverse(node, depth=0, is_main_pathway=True):
        nonlocal max_depth, heterocycle_molecules, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and node.get("smiles"):
            # Skip starting materials
            if node.get("in_stock", False):
                return

            # Only track molecules in the main synthetic pathway
            if is_main_pathway:
                main_pathway_molecules.append(node["smiles"])

                # Check for various heterocycles
                for ring in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, node["smiles"]):
                        heterocycle_molecules.append((node["smiles"], ring))
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        break

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical' (or 'mol'), depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth, is_main_pathway)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a valid heterocycle preservation strategy:
    # 1. At least 2 molecules in the main pathway have heterocycles
    # 2. The synthesis has some depth (not just a single reaction)
    if len(heterocycle_molecules) >= 2 and max_depth >= 1:
        print(
            f"Detected heterocycle preservation strategy with {len(heterocycle_molecules)} heterocycle-containing molecules"
        )
        print(f"Heterocycles found: {heterocycle_molecules}")

        # Record structural constraint: at least 2 molecules with target heterocycle
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "molecule_with_target_heterocycle",
                "operator": ">=",
                "value": 2
            }
        })
        # Record structural constraint: at least 1 reaction step (max_depth >= 1 implies this)
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_steps",
                "operator": ">=",
                "value": 1
            }
        })

        # Check if the target molecule (root) has a heterocycle
        target_has_heterocycle = False
        if route["type"] == "mol" and route.get("smiles"):
            for ring in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring, route["smiles"]):
                    target_has_heterocycle = True
                    break

        # For a true heterocycle preservation strategy, the target molecule must have a heterocycle
        if target_has_heterocycle or (
            heterocycle_molecules and route["smiles"] == heterocycle_molecules[0][0]
        ):
            result = True
            # Record structural constraint: final product has target heterocycle
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "any_target_heterocycle",
                    "position": "final_product"
                }
            })

    return result, findings_json
