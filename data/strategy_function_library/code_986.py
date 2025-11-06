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
    This function detects if the synthesis involves a beta-lactam core structure
    that is preserved throughout the synthesis.
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

    # Beta-lactam is a 4-membered ring with N and C=O
    # Define the beta-lactam pattern correctly
    beta_lactam_pattern = Chem.MolFromSmarts("[#6]1[#6](=[O])[#7][#6]1")

    # Track if we've found a beta-lactam that's preserved
    beta_lactam_preserved = [False]

    def check_beta_lactam(smiles):
        """Helper function to check if a molecule contains a beta-lactam"""
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(beta_lactam_pattern):
            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal beta_lactam_preserved, findings_json
        # Check molecule nodes for beta-lactam
        if node["type"] == "mol" and "smiles" in node:
            if check_beta_lactam(node["smiles"]):
                print(f"Beta-lactam found in molecule at depth {depth}: {node['smiles']}")
                findings_json["atomic_checks"]["ring_systems"].append("beta-lactam")

                # If this is the final product (depth 0), mark as preserved initially
                if depth == 0:
                    beta_lactam_preserved[0] = True
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "beta-lactam",
                            "position": "last_stage"
                        }
                    })

        # Check reaction nodes to see if beta-lactam is preserved
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_has_beta_lactam = check_beta_lactam(product)
                reactant_has_beta_lactam = any(check_beta_lactam(r) for r in reactants)

                # If a reactant has a beta-lactam but the product doesn't, it's destroyed
                if not product_has_beta_lactam and reactant_has_beta_lactam:
                    print(f"Beta-lactam destroyed in reaction at depth {depth}: {rsmi}")
                    beta_lactam_preserved[0] = False
                    # If beta-lactam is destroyed, it means the negation constraint is NOT met
                    # We only add to findings_json if the constraint IS met.
                    # So, if beta_lactam_preserved[0] becomes False, we don't add the negation constraint.
                elif product_has_beta_lactam and reactant_has_beta_lactam:
                    # If beta-lactam is present in both reactant and product, it's preserved through this step
                    pass # No specific finding to add here, preservation is cumulative
                elif not reactant_has_beta_lactam and product_has_beta_lactam:
                    # Beta-lactam formed in this step, but not necessarily preserved from earlier
                    pass

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Continue traversing children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not 'reaction'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    result = beta_lactam_preserved[0]

    # Add the negation constraint if beta-lactam was preserved throughout
    if result:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_destruction"
            }
        })

    return result, findings_json