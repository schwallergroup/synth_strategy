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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Ester with primary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving late-stage amide coupling.
    It specifically checks if the final reaction step (depth=1) is one of the named reactions
    defined in the AMIDE_COUPLING_REACTIONS list.
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

    amide_coupling_detected = False
    final_product_smiles = route["smiles"]

    print(f"Analyzing route for late-stage amide coupling. Final product: {final_product_smiles}")

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_detected, findings_json

        # For reaction nodes, check if it's an amide coupling reaction
        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi") and depth <= 1:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            parts = rsmi.split(">")
            if len(parts) >= 3:  # Ensure proper format
                # Check if this is an amide coupling reaction
                is_amide_coupling = False

                # Try different amide coupling reaction types
                for reaction_type in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected amide coupling reaction: {reaction_type}")
                        is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                # If we found amide coupling at a late stage (depth <= 1)
                if is_amide_coupling:
                    print(f"Late-stage amide coupling confirmed at depth {depth}")
                    amide_coupling_detected = True
                    # Add the structural constraint if detected
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": [
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                "Carboxylic acid with primary amine to amide",
                                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                "Ester with primary amine to amide",
                                "Acylation of primary amines",
                                "Acylation of secondary amines",
                                "Schotten-Baumann_amide"
                            ],
                            "position": "last_stage"
                        }
                    })

        # Continue traversing
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: amide_coupling_detected = {amide_coupling_detected}")
    return amide_coupling_detected, findings_json
