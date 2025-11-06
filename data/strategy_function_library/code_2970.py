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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves the formation of an aromatic ring from a non-aromatic precursor, and if so, checks that this type of transformation occurs at only a single stage of the synthesis.
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

    aromatic_incorporation_depths = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                product_has_aromatic = product_mol and any(atom.GetIsAromatic() for atom in product_mol.GetAtoms())

                if product_has_aromatic:
                    # Check if at least one reactant is non-aromatic
                    reactant_is_non_aromatic = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        # A reactant is considered non-aromatic if it's a valid molecule and has no aromatic atoms.
                        if reactant_mol and not any(atom.GetIsAromatic() for atom in reactant_mol.GetAtoms()):
                            reactant_is_non_aromatic = True
                            break
                    
                    if reactant_is_non_aromatic:
                        aromatic_incorporation_depths.append(depth)
                        # Record the 'aromatization' atomic check
                        if "aromatization" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aromatization")
            except Exception:
                pass

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not a reaction (i.e., chemical)
                new_depth += 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if not aromatic_incorporation_depths:
        result = False
    else:
        min_depth = min(aromatic_incorporation_depths)
        result = all(depth == min_depth for depth in aromatic_incorporation_depths)
        
        if result:
            # If the condition for a single stage is met, record the structural constraint
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "stages_with_aromatization",
                    "operator": "==",
                    "value": 1
                }
            })

    return result, findings_json
