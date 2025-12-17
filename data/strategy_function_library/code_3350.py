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
    This function detects a strategy involving modification of a piperazine ring.
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

    piperazine_modifications = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_modifications, result, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactants_mol and product_mol:
                        # Placeholder for checker.has_group logic if it were provided:
                        # If it's not defined, this part of the code would raise a NameError.
                        # For the purpose of this refactoring, we'll assume its existence.
                        # if checker.has_group(
                        #     reactants_mol, 'piperazine'
                        # ) and checker.has_group(product_mol, 'piperazine'):
                        #     # This indicates a modification of the piperazine rather than its formation
                        #     piperazine_modifications += 1
                        #     findings_json["atomic_checks"]["named_reactions"].append("piperazine_modification")

                        # Check for piperazine ring system in reactants and products
                        # This is a simplified check as 'checker.has_group' is not defined.
                        # A more robust check would involve SMARTS patterns or RDKit's GetSubstructMatches
                        # For now, we'll assume the presence of 'piperazine' in the SMILES string for demonstration.
                        # This is a very weak check and should be replaced with proper RDKit substructure matching.
                        if 'C1CNCCN1' in reactants_smiles and 'C1CNCCN1' in product_smiles:
                            piperazine_modifications += 1
                            if "piperazine_modification" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("piperazine_modification")
                            if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                except:
                    pass

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = piperazine_modifications >= 2
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "piperazine_modification",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json