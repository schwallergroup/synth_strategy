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


CN_COUPLING_REACTIONS = [
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Goldberg coupling",
    "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride",
    "Ullmann-Goldberg Substitution amine",
    "N-arylation_heterocycles",
    "Buchwald-Hartwig",
    "Chan-Lam amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final or penultimate step) aromatic C-N bond formation
    by checking if the reaction is one of the following types: Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine,
    Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine, N-arylation (Buchwald-Hartwig/Ullmann-Goldberg),
    Goldberg coupling, Goldberg coupling aryl amine-aryl chloride, Goldberg coupling aryl amide-aryl chloride,
    Ullmann-Goldberg Substitution amine, N-arylation_heterocycles, Buchwald-Hartwig, or Chan-Lam amine.
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

    found_late_cn_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_cn_formation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            node_depth = depth

            # Check late-stage reactions (depth 1 is the final step).
            if 1 <= node_depth <= 2:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a known C-N bond forming reaction.
                for reaction_name in CN_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        found_late_cn_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        # If a reaction is found, and it's in the correct position, add the structural constraint
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": [
                                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                                    "Goldberg coupling",
                                    "Goldberg coupling aryl amine-aryl chloride",
                                    "Goldberg coupling aryl amide-aryl chloride",
                                    "Ullmann-Goldberg Substitution amine",
                                    "N-arylation_heterocycles",
                                    "Buchwald-Hartwig",
                                    "Chan-Lam amine"
                                ],
                                "position": "final_or_penultimate_step"
                            }
                        })
                        break # Break from the inner loop once a reaction is found

        # Continue traversing.
        for child in node.get("children", []):
            # Optimization: stop traversing if we already found the flag.
            if found_late_cn_formation:
                break

            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal.
    dfs_traverse(route)
    return found_late_cn_formation, findings_json
