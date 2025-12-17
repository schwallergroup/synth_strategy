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
    "isoxazole",
    "oxazole",
    "thiazole",
    "pyrrole",
    "furan",
    "pyridine",
    "imidazole",
    "triazole",
    "tetrazole",
    "pyrazole",
    "oxadiazole",
    "thiadiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a heterocycle from a predefined list is formed in the middle stages of a synthesis.
    The middle stage is defined as occurring between 25% and 75% of the total synthesis depth.
    The list of heterocycles checked includes: isoxazole, oxazole, thiazole, pyrrole, furan, pyridine, imidazole, triazole, tetrazole, pyrazole, oxadiazole, thiadiazole, pyrimidine, pyrazine, and pyridazine.
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

    # Track if we found heterocycle formation
    heterocycle_formed = False
    total_depth = 0
    max_depth = 0
    heterocycle_depth = None

    # List of heterocycles to check for
    heterocycles = HETEROCYCLES_OF_INTEREST

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed, total_depth, max_depth, heterocycle_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check which heterocycles are formed in this reaction
            formed_heterocycles = []
            for ring in heterocycles:
                # Check if ring is in product but not in any reactant
                if checker.check_ring(ring, product_smiles):
                    if not any(checker.check_ring(ring, r) for r in reactants_smiles):
                        formed_heterocycles.append(ring)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

            # If any heterocycle is formed in this step
            if formed_heterocycles:
                heterocycle_formed = True
                heterocycle_depth = depth
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    total_depth = max_depth

    result = False
    # Check if heterocycle formation occurred in the middle stages
    # (not at the beginning or end of the synthesis)
    if heterocycle_formed and heterocycle_depth is not None:
        # Consider middle stage as between 25% and 75% of the total synthesis depth
        if total_depth <= 1: # Avoid division by zero and handle short syntheses
            result = False
        
        lower_bound = total_depth * 0.25
        upper_bound = total_depth * 0.75

        is_mid_stage = lower_bound <= heterocycle_depth <= upper_bound

        if is_mid_stage:
            result = True
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "ring_formation",
                    "position": "middle_stage (25%-75% of total depth)"
                }
            })

    return result, findings_json