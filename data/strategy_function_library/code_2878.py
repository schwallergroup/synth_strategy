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


ESTER_HYDROLYSIS_REACTIONS = [
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
]

ESTERIFICATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "Protection of carboxylic acid",
    "O-alkylation of carboxylic acids with diazo compounds",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a bidirectional functional group interconversion strategy involving an ester and a carboxylic acid, such as an ester → acid → ester or acid → ester → acid sequence. The function identifies individual steps by checking against predefined lists of named reactions for ester hydrolysis (ESTER_HYDROLYSIS_REACTIONS) and esterification (ESTERIFICATION_REACTIONS).
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
    result = False

    # Track transformations in sequence with molecule information
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product_smiles = rsmi.split(">")[-1]

                # Check for ester hydrolysis (ester → acid)
                for r in ESTER_HYDROLYSIS_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        transformations.append(("ester_to_acid", product_smiles, depth))
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break # Only add one if found

                # Check for esterification (acid → ester)
                for r in ESTERIFICATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        transformations.append(("acid_to_ester", product_smiles, depth))
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break # Only add one if found

            except Exception:
                # Silently ignore errors in reaction processing to not halt the entire analysis.
                pass

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Sort transformations by depth to ensure correct sequence analysis
    transformations.sort(key=lambda x: x[2])

    # Check for the sequence ester→acid→ester in the synthetic route
    for i in range(len(transformations) - 1):
        if (
            transformations[i][0] == "ester_to_acid"
            and transformations[i + 1][0] == "acid_to_ester"
        ):
            result = True
            # Add structural constraint for ester->acid->ester sequence
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "ester_to_acid",
                        "acid_to_ester"
                    ],
                    "is_contiguous": True
                }
            })
            break # Found one, no need to check further for this sequence

    # Also check for the reverse sequence (acid→ester→acid)
    # Note: If both sequences are found, 'result' will be True, and both structural constraints will be added.
    for i in range(len(transformations) - 1):
        if (
            transformations[i][0] == "acid_to_ester"
            and transformations[i + 1][0] == "ester_to_acid"
        ):
            result = True
            # Add structural constraint for acid->ester->acid sequence
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "acid_to_ester",
                        "ester_to_acid"
                    ],
                    "is_contiguous": True
                }
            })
            break # Found one, no need to check further for this sequence

    return result, findings_json
