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
    Detects if the synthetic route contains multiple Boc deprotection steps.
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

    boc_deprotection_count = 0
    boc_pattern = Chem.MolFromSmarts("[C][C]([C])([C])[O][C](=[O])[N]")

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_count, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                products_mol = Chem.MolFromSmiles(products_smiles)

                if reactants_mol and products_mol:
                    # Check if Boc group is in reactants but not in products
                    if reactants_mol.HasSubstructMatch(
                        boc_pattern
                    ):
                        findings_json["atomic_checks"]["functional_groups"].append("Boc group")
                        if not products_mol.HasSubstructMatch(boc_pattern):
                            boc_deprotection_count += 1
                            findings_json["atomic_checks"]["named_reactions"].append("Boc deprotection")
            except:
                # Error processing is silenced for robustness in a large-scale pipeline
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    result = boc_deprotection_count >= 2

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "Boc deprotection",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
