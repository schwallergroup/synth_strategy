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


HETEROCYCLES_OF_INTEREST = [
    "imidazole",
    "triazole",
    "tetrazole",
    "oxazole",
    "thiazole",
    "pyrazole",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
    "pyrrole",
    "furan",
    "thiophene",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "indole",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis route involves the late-stage formation of a new heterocyclic ring.
    The function checks for the de novo formation of any heterocycle specified in the HETEROCYCLES_OF_INTEREST list.
    A formation event is considered 'late-stage' if it occurs in the final half of the total synthesis steps (e.g., in steps 1-5 of a 10-step synthesis).
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

    heterocycle_formation_depth = None
    max_depth = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        # Only check for formation if we haven't found an event yet.
        # We must continue traversing regardless to find the true max_depth.
        if heterocycle_formation_depth is None:
            if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Identify all heterocycles of interest in the product
                product_heterocycles = [
                    h for h in HETEROCYCLES_OF_INTEREST if checker.check_ring(h, product_str)
                ]

                # If any are found, check they are not present in the reactants
                if product_heterocycles:
                    reactant_heterocycles = {
                        h
                        for reactant in reactants_str.split(".")
                        for h in HETEROCYCLES_OF_INTEREST
                        if checker.check_ring(h, reactant)
                    }

                    # A new heterocycle is formed if any in the product are not in the reactants
                    newly_formed_heterocycles = [h for h in product_heterocycles if h not in reactant_heterocycles]
                    if newly_formed_heterocycles:
                        heterocycle_formation_depth = depth
                        # Record atomic checks for newly formed heterocycles
                        for h in newly_formed_heterocycles:
                            if h not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(h)
                        # Record the conceptual 'ring_formation' named reaction
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Recursively traverse all children to find the full depth of the synthesis
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if heterocycle formation occurs in the late stage (first half of synthesis)
    if heterocycle_formation_depth is not None and max_depth > 0:
        # Lower depth values correspond to later stages in synthesis
        is_late_stage = heterocycle_formation_depth <= (max_depth / 2)
        if is_late_stage:
            result = True
            # Add structural constraint if late-stage formation is detected
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "ring_formation",
                    "position": "final_half"
                }
            })

    return result, findings_json
