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
    """
    This function detects if the synthesis involves the convergent assembly
    of three key fragments: methylphenyl core, piperidine scaffold, and cyano-pyridine.
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

    # Track fragments and their positions in the synthesis
    fragment_paths = {"methylphenyl": [], "piperidine": [], "cyanopyridine": []}

    # Track depth of each node to identify convergent steps
    node_depths = {}

    def identify_fragments(mol_smiles):
        """Identify which fragments are present in a molecule"""
        fragments = set()

        # Check for methylphenyl core (phenyl with methyl or connected to triazole/imidazole)
        if checker.check_ring("benzene", mol_smiles):
            findings_json["atomic_checks"]["ring_systems"].append("benzene")
            if any(
                [
                    checker.check_ring("triazole", mol_smiles),
                    checker.check_ring("imidazole", mol_smiles),
                    checker.check_fg("Methyl", mol_smiles),
                ]
            ):
                if checker.check_ring("triazole", mol_smiles):
                    findings_json["atomic_checks"]["ring_systems"].append("triazole")
                if checker.check_ring("imidazole", mol_smiles):
                    findings_json["atomic_checks"]["ring_systems"].append("imidazole")
                if checker.check_fg("Methyl", mol_smiles):
                    findings_json["atomic_checks"]["functional_groups"].append("Methyl")
                fragments.add("methylphenyl")

        # Check for piperidine scaffold
        if checker.check_ring("piperidine", mol_smiles):
            findings_json["atomic_checks"]["ring_systems"].append("piperidine")
            fragments.add("piperidine")

        # Check for cyano-pyridine
        if checker.check_ring("pyridine", mol_smiles):
            findings_json["atomic_checks"]["ring_systems"].append("pyridine")
            if checker.check_fg("Nitrile", mol_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                fragments.add("cyanopyridine")

        return fragments

    def dfs_traverse(node, depth=0):
        # Store the depth of this node
        node_id = id(node)
        node_depths[node_id] = depth

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            fragments = identify_fragments(mol_smiles)

            # Record which fragments are found at this node
            for fragment in fragments:
                fragment_paths[fragment].append((node_id, depth))

            # Check if this is a starting material
            if node.get("in_stock", False):
                print(f"Found starting material with fragments: {fragments}, SMILES: {mol_smiles}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # node['type'] == 'mol'
                dfs_traverse(child, depth + 1)

    # Traverse the route to identify fragments
    dfs_traverse(route)

    # Count how many unique fragments were found
    fragments_found = sum(1 for paths in fragment_paths.values() if paths)
    print(f"Found {fragments_found} unique fragments")

    # Check if we have all three fragments
    if fragments_found < 3:
        # Try to identify the missing fragment from the starting materials
        missing_fragments = [f for f, paths in fragment_paths.items() if not paths]
        print(f"Missing fragments: {missing_fragments}")
        return False, findings_json
    else:
        # Record the structural constraint if 3 or more fragments are found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_fragments_found",
                "operator": ">=",
                "value": 3,
                "description": "Counts three specific fragments: 1) methylphenyl (benzene + (triazole|imidazole|methyl)), 2) piperidine, 3) cyanopyridine (pyridine + nitrile)."
            }
        })

    # Analyze if the synthesis is convergent by checking depths
    # For convergent synthesis, fragments should be introduced early (higher depth)
    # and combined later (lower depth)

    # Get the earliest (highest depth) appearance of each fragment
    earliest_appearances = {}
    for fragment, paths in fragment_paths.items():
        if paths:
            # Sort by depth in descending order (highest depth first)
            sorted_paths = sorted(paths, key=lambda x: x[1], reverse=True)
            earliest_appearances[fragment] = sorted_paths[0]

    # Check if fragments appear at different branches (different node IDs)
    unique_nodes = len(set(node_id for node_id, _ in earliest_appearances.values()))

    print(f"Fragments appear at {unique_nodes} different branches")
    print(f"Earliest appearances: {earliest_appearances}")

    # For convergent synthesis:
    # 1. We need all three fragments
    # 2. They should appear at different branches (at least 2 different branches)
    # 3. They should be combined later in the synthesis

    is_convergent = fragments_found >= 3 and unique_nodes >= 2

    if unique_nodes >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "distinct_synthesis_branches_for_fragments",
                "operator": ">=",
                "value": 2,
                "description": "The earliest appearances of the required fragments must originate from at least 2 distinct synthesis branches."
            }
        })

    if is_convergent:
        print("Found convergent assembly of three fragments")

    return is_convergent, findings_json
