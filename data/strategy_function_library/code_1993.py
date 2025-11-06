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


TETRAMETHYLTETRALIN_ISOMER_SMARTS = [
    Chem.MolFromSmarts("C1CCc2c(C)c(C)c(C)c(C)c2C1"),
    Chem.MolFromSmarts("CC1(C)CCC(C)(C)c2ccccc21"),
    Chem.MolFromSmarts("C1CC(C)(C)c2ccc(C)c(C)c2C1"),
    Chem.MolFromSmarts("C1CC(C)(C)c2cc(C)c(C)cc2C1"),
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final product of a synthesis contains a specific tetramethyltetralin isomer.
    The isomers are defined by the SMARTS patterns in the module-level list
    TETRAMETHYLTETRALIN_ISOMER_SMARTS.
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

    scaffold_preserved = True
    found_in_final_product = False

    def dfs_traverse(node, depth=0):
        nonlocal scaffold_preserved, found_in_final_product, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol:
                # If it's the final product (depth 0), we must find the scaffold
                if depth == 0:
                    scaffold_found = False
                    for pattern in TETRAMETHYLTETRALIN_ISOMER_SMARTS:
                        if mol.HasSubstructMatch(pattern):
                            scaffold_found = True
                            found_in_final_product = True
                            findings_json["atomic_checks"]["ring_systems"].append("tetramethyltetralin isomer")
                            break

                    if not scaffold_found:
                        scaffold_preserved = False

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node.get("type") != "reaction": # If current node is not 'reaction' (e.g., 'chemical' or 'mol')
            next_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    dfs_traverse(route)

    if not found_in_final_product:
        scaffold_preserved = False
    else:
        # If found_in_final_product is True, it means the scaffold was found in the last stage.
        # Add the structural constraint.
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "tetramethyltetralin isomer",
                "position": "last_stage"
            }
        })

    return scaffold_preserved, findings_json
