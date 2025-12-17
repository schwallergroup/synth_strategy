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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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
    Detects preservation of difluoromethoxy group throughout the synthesis
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

    # Track if difluoromethoxy group is present in the main synthetic pathway
    all_steps_have_difluoromethoxy = True
    step_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal all_steps_have_difluoromethoxy, step_count, findings_json

        if node["type"] == "mol" and "smiles" in node:
            # Only check molecules that are part of the main synthetic pathway
            # (not in_stock reagents or catalysts)
            if depth == 0 or (node.get("children") and len(node.get("children")) > 0):
                mol_smiles = node["smiles"]

                # Check for difluoromethoxy (OCF2H)
                difluoromethoxy_pattern = Chem.MolFromSmarts("OC(F)F")
                mol = Chem.MolFromSmiles(mol_smiles)
                has_difluoromethoxy = mol and mol.HasSubstructMatch(difluoromethoxy_pattern)

                if has_difluoromethoxy:
                    # If difluoromethoxy is found, record it
                    if "difluoromethoxy" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("difluoromethoxy")
                else:
                    all_steps_have_difluoromethoxy = False

                step_count += 1

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or 'chemical' type for other cases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Only return True if we've checked at least one molecule and all have difluoromethoxy
    result = all_steps_have_difluoromethoxy and step_count > 0

    # Populate structural constraints based on the final result
    if not all_steps_have_difluoromethoxy:
        # This corresponds to the negation constraint: "absence_of_difluoromethoxy"
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "absence_of_difluoromethoxy",
                "scope": "main_pathway_molecules"
            }
        })
    
    if step_count > 0:
        # This corresponds to the count constraint: "main_pathway_molecule" > 0
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "main_pathway_molecule",
                "operator": ">",
                "value": 0
            }
        })

    return result, findings_json
