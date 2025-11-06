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


N_ALKYLATION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Alkylation of amines",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves sequential C-N bond formations:
    first an N-alkylation, then an amide bond formation.
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
    found_amide_formation = False
    n_alkylation_depth = -1
    amide_formation_depth = -1

    print("Starting sequential C-N bond formation strategy check")

    def dfs_traverse(node, current_depth=0):
        nonlocal found_n_alkylation, found_amide_formation, n_alkylation_depth, amide_formation_depth, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Extract depth from metadata if available
                depth = current_depth
                if "Depth" in node.get("metadata", {}):
                    try:
                        depth_value = node["metadata"]["Depth"]
                        if isinstance(depth_value, (int, float)):
                            depth = int(depth_value)
                        elif isinstance(depth_value, str) and depth_value.isdigit():
                            depth = int(depth_value)
                    except (ValueError, TypeError):
                        pass

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for N-alkylation reactions
                for rxn in N_ALKYLATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        found_n_alkylation = True
                        n_alkylation_depth = depth
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Found N-alkylation at depth {depth}")
                        break # Found one, no need to check others for this reaction node

                # Check for amide formation reactions
                for rxn in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        found_amide_formation = True
                        amide_formation_depth = depth
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Found amide formation at depth {depth}")
                        break # Found one, no need to check others for this reaction node

        # Traverse children with modified depth
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, current_depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both reactions were found and in the correct order
    # Remember: higher depth = earlier in synthesis (retrosynthetic direction)
    result = (
        found_n_alkylation and found_amide_formation and n_alkylation_depth > amide_formation_depth
    )
    print(f"N-alkylation found: {found_n_alkylation} at depth {n_alkylation_depth}")
    print(f"Amide formation found: {found_amide_formation} at depth {amide_formation_depth}")
    print(f"Correct sequence (N-alkylation before amide formation): {result}")

    if found_n_alkylation and found_amide_formation:
        # Add co-occurrence constraint if both are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "N-alkylation",
                    "amide_formation"
                ]
            }
        })

    if result:
        # Add sequence constraint if the order is correct
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "N-alkylation",
                "after": "amide_formation"
            }
        })

    return result, findings_json
