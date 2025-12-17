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
    Detects if a thiophene or thiazole heterocycle is present in every molecule of the main synthetic pathway.
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

    # Track if the main synthetic pathway preserves thiophene or thiazole
    thiophene_in_main_pathway = True

    # Track visited molecules to avoid checking reagents multiple times
    visited_molecules = set()

    def dfs_traverse(node, depth=0, is_main_path=True):
        nonlocal thiophene_in_main_pathway, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Skip if we've already processed this molecule
            if mol_smiles in visited_molecules:
                return

            visited_molecules.add(mol_smiles)

            # Only check molecules in the main synthetic pathway
            if is_main_path:
                # Check if the molecule contains a thiophene or thiazole ring
                has_thiophene = checker.check_ring("thiophene", mol_smiles)
                has_thiazole = checker.check_ring("thiazole", mol_smiles)

                if has_thiophene:
                    if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("thiophene")
                if has_thiazole:
                    if "thiazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("thiazole")

                if not (has_thiophene or has_thiazole):
                    thiophene_in_main_pathway = False

        # Determine the depth for the next recursive call based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Depth increases when moving from chemical to reaction
            next_depth = depth + 1

        # For reaction nodes, determine which children are part of the main pathway
        if node["type"] == "reaction" and "children" in node:
            # In retrosynthetic analysis, the first child is typically the product
            # and the rest are reactants/reagents
            for i, child in enumerate(node.get("children", [])):
                # In retrosynthesis, the first child (i==0) is the product
                # and is part of the main pathway
                child_is_main_path = is_main_path and i == 0
                dfs_traverse(child, next_depth, child_is_main_path)

        # Traverse other children (e.g., from mol nodes)
        else:
            for child in node.get("children", []):
                dfs_traverse(child, next_depth, is_main_path)

    # Start traversal from the root
    dfs_traverse(route)

    if thiophene_in_main_pathway:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "main_path_molecule_without_thiophene_and_thiazole",
                "operator": "==",
                "value": 0
            }
        })

    return thiophene_in_main_pathway, findings_json
