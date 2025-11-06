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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves indazole ring formation.
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

    indazole_formation_detected = False

    def dfs_traverse(node, depth):
        nonlocal indazole_formation_detected, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if indazole is in product and not in any reactants
                if checker.check_ring("indazole", product_smiles):
                    findings_json["atomic_checks"]["ring_systems"].append("indazole")
                    reactants_have_indazole = False
                    for reactant_smiles in reactants_smiles:
                        if checker.check_ring("indazole", reactant_smiles):
                            reactants_have_indazole = True
                            break

                    if not reactants_have_indazole:
                        indazole_formation_detected = True
                        # The strategy JSON implies 'ring_formation' as a named reaction
                        # This is a conceptual mapping as there's no explicit 'check_reaction' for it.
                        # Assuming this condition implies 'ring_formation'.
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
            except Exception:
                # Silently ignore errors in single reaction processing
                pass

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Chemical node to reaction node
                new_depth = depth + 1
            # If node['type'] is 'reaction', new_depth remains 'depth' (reaction node to chemical node)
            dfs_traverse(child, new_depth)

    def traverse(node):
        # The recursive traversal logic
        dfs_traverse(node, 0) # Start initial depth at 0

    traverse(route)
    return indazole_formation_detected, findings_json