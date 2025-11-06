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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route includes a multi-component reaction
    (a reaction where 3 or more reactant fragments combine to form fewer product fragments).
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

    mcr_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal mcr_detected, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Count number of disconnected fragments in reactants
                reactant_fragments = reactants_part.split(".")

                # Count number of disconnected fragments in product
                product_fragments = product_part.split(".")

                # If 3+ reactant fragments combine to form fewer product fragments
                if len(reactant_fragments) >= 3 and len(product_fragments) < len(
                    reactant_fragments
                ):
                    mcr_detected = True
                    findings_json["atomic_checks"]["named_reactions"].append("multi-component_reaction")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return mcr_detected, findings_json