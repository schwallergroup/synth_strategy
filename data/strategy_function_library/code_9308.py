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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a piperazine moiety is
    incorporated in a late-stage reaction.
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

    # Track if we found the pattern
    found_piperazine_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_piperazine_incorporation, findings_json

        # Add depth to metadata if not present
        if node["type"] == "reaction":
            if "metadata" not in node:
                node["metadata"] = {}
            if "depth" not in node["metadata"]:
                node["metadata"]["depth"] = depth

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "depth" in node["metadata"]
            and node["metadata"]["depth"] <= 1
        ):
            # This is a late-stage reaction (depth 0 or 1)
            # Add structural constraint for depth
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "ring_formation",
                    "position": "depth <= 1"
                }
            })
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has piperazine
                has_piperazine_product = checker.check_ring("piperazine", product_smiles)
                if has_piperazine_product:
                    findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                if has_piperazine_product:
                    # Check if any reactant has piperazine
                    reactants_with_piperazine = [
                        checker.check_ring("piperazine", r) for r in reactants_smiles if r
                    ]

                    # True incorporation means piperazine is formed, not carried through
                    if not any(reactants_with_piperazine):
                        found_piperazine_incorporation = True
                        # Inferring 'ring_formation' from the logic where piperazine is in product but not reactants
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
            except Exception as e:
                pass

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_piperazine_incorporation, findings_json
