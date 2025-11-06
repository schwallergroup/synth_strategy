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
    Detects the formation of a biphenyl scaffold.
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

    # Track if we found biphenyl construction
    found_biphenyl_construction = False

    def is_biphenyl(smiles):
        """Check if a molecule contains a biphenyl scaffold"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")
            if pattern and mol.HasSubstructMatch(pattern):
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_biphenyl_construction, findings_json

        if node["type"] == "reaction":
            reaction_smiles = node["metadata"].get("rsmi", "")
            if reaction_smiles:
                # Check if product contains biphenyl but reactants don't
                reactants = reaction_smiles.split(">")[0].split(".")
                product = reaction_smiles.split(">")[-1]

                if is_biphenyl(product):
                    findings_json["atomic_checks"]["ring_systems"].append("biphenyl")
                    # Check if any reactant already has the biphenyl scaffold
                    reactant_has_biphenyl = any(is_biphenyl(r) for r in reactants)

                    if not reactant_has_biphenyl:
                        found_biphenyl_construction = True
                        findings_json["atomic_checks"]["named_reactions"].append("biphenyl_formation")
                        print(f"Found biphenyl construction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if found_biphenyl_construction:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "biphenyl_formation",
                "operator": ">=",
                "value": 1
            }
        })

    print(f"Biphenyl scaffold construction strategy: {found_biphenyl_construction}")
    return found_biphenyl_construction, findings_json
