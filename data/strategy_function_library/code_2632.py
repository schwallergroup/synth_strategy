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


N_ALKYLATION_REACTION_TYPES = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "DMS Amine methylation",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
    "N-methylation",
    "Alkylation of amines",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage N-alkylation in the final step of the synthesis by checking for a specific set of named reactions.
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

    found_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for N-alkylation reactions using the checker
                for reaction_type in N_ALKYLATION_REACTION_TYPES:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected late-stage N-alkylation: {reaction_type} at depth {depth}")
                        found_n_alkylation = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if found at the correct depth
                        if depth == 0: # Assuming 'last_stage' means depth 0 (the final reaction)
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "targets": [
                                        "N-alkylation of primary amines with alkyl halides",
                                        "N-alkylation of secondary amines with alkyl halides",
                                        "Methylation with MeI_primary",
                                        "Methylation with MeI_secondary",
                                        "Methylation with MeI_tertiary",
                                        "DMS Amine methylation",
                                        "Eschweiler-Clarke Primary Amine Methylation",
                                        "Eschweiler-Clarke Secondary Amine Methylation",
                                        "Reductive methylation of primary amine with formaldehyde",
                                        "N-methylation",
                                        "Alkylation of amines",
                                        "Reductive amination with aldehyde",
                                        "Reductive amination with ketone",
                                        "Reductive amination with alcohol"
                                    ],
                                    "position": "last_stage"
                                }
                            })
                        return  # Exit early once found

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

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

    print(f"Late-stage N-alkylation detected: {found_n_alkylation}")
    return found_n_alkylation, findings_json
