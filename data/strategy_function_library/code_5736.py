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


THIOETHER_COUPLING_REACTIONS = [
    "S-alkylation of thiols",
    "S-alkylation of thiols with alcohols",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where a pyridine ring is formed, followed by a late-stage fragment coupling via a specific thioether-forming reaction. The specific coupling reactions checked are 'S-alkylation of thiols' and 'S-alkylation of thiols with alcohols'.
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

    # Track if we found the key features and their depths
    heterocycle_formation_depth = None
    thioether_coupling_depth = None

    result = False # Initialize the boolean result

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, thioether_coupling_depth, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for heterocycle formation (pyridine ring)
                # First check if product contains pyridine
                if checker.check_ring("pyridine", product_smiles):
                    # Check if any reactant doesn't have pyridine
                    if not all(checker.check_ring("pyridine", r) for r in reactants_smiles):
                        print(f"Found heterocycle formation (pyridine ring) at depth {depth}")
                        # Update depth only if not already set or if this is earlier (higher depth)
                        if (
                            heterocycle_formation_depth is None
                            or depth > heterocycle_formation_depth
                        ):
                            heterocycle_formation_depth = depth
                            if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # Check for thioether coupling (S-alkylation)
                # Check if this is a specific S-alkylation reaction
                for rxn in THIOETHER_COUPLING_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Found thioether coupling via S-alkylation at depth {depth}")
                        # Update depth only if not already set or if this is later (lower depth)
                        if thioether_coupling_depth is None or depth < thioether_coupling_depth:
                            thioether_coupling_depth = depth
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break # Found one, no need to check other reactions

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if both key features were found and thioether coupling is at a later stage (lower depth)
    if heterocycle_formation_depth is not None and thioether_coupling_depth is not None:
        print(
            f"Heterocycle formation depth: {heterocycle_formation_depth}, Thioether coupling depth: {thioether_coupling_depth}"
        )
        if heterocycle_formation_depth > thioether_coupling_depth:
            result = True
            # Add the structural constraint if the condition is met
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "ring_formation",
                    "after": [
                        "S-alkylation of thiols",
                        "S-alkylation of thiols with alcohols"
                    ]
                }
            })

    return result, findings_json
