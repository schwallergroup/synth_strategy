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


N_ALKYLATION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "Methylation with MeI_aryl",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage N-dealkylation strategy by identifying a specific N-alkylation reaction in the final synthetic step.
    The forward-direction reaction is checked against a defined list of types, including various methylations (e.g., with MeI, Eschweiler-Clarke) and general N-alkylations with alkyl halides.
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

    found_n_dealkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_dealkylation, findings_json
        if found_n_dealkylation:
            return

        # Per problem spec, depth=1 is the final step.
        is_late_stage = depth <= 1

        if node.get("type") == "reaction" and is_late_stage and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # In the forward direction, we're looking for N-alkylation reactions,
            # which correspond to N-dealkylation in retrosynthesis.
            for rxn in N_ALKYLATION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    found_n_dealkylation = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    # Add the structural constraint when the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": [
                                "N-alkylation of primary amines with alkyl halides",
                                "N-alkylation of secondary amines with alkyl halides",
                                "Methylation with MeI_primary",
                                "Methylation with MeI_secondary",
                                "Methylation with MeI_tertiary",
                                "Methylation with MeI_aryl",
                                "Eschweiler-Clarke Primary Amine Methylation",
                                "Eschweiler-Clarke Secondary Amine Methylation"
                            ],
                            "position": "last_stage"
                        }
                    })
                    return

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node.get("type") == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_n_dealkylation, findings_json
