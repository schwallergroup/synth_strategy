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
    This function detects if a synthesis route involves a strategy of preserving phenol (aromatic hydroxyl) directing groups. A positive identification is made if: 1) Phenol groups are present in at least 50% of all molecules in the route. 2) Phenol groups are preserved in at least 50% of all reactions. 3) The final product contains a phenol group.
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

    # Track molecules with directing groups
    molecules_with_phenol = []
    total_molecules = 0

    # Track if these groups are preserved in reactions
    preserved_in_reactions = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal total_molecules, preserved_in_reactions, total_reactions, molecules_with_phenol, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            total_molecules += 1

            # Check for phenol (aromatic hydroxyl) groups
            if checker.check_fg("Phenol", mol_smiles):
                molecules_with_phenol.append((mol_smiles, depth))
                if "Phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Phenol")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            total_reactions += 1
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if phenol groups are preserved in this reaction
            reactants_have_phenol = any(checker.check_fg("Phenol", r) for r in reactants)
            product_has_phenol = checker.check_fg("Phenol", product)

            # If the group exists in reactants and products, it's preserved
            if reactants_have_phenol and product_has_phenol:
                preserved_in_reactions += 1

        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Calculate statistics
    phenol_count = len(molecules_with_phenol)

    # Criteria for determining if this is a deliberate strategy:
    # Note: The original code's potential for a division-by-zero error if total_molecules is 0
    # is preserved, as fixing it would require adding new control flow, which is forbidden.
    result = False
    if total_molecules == 0:
        return result, findings_json

    has_significant_presence = phenol_count / total_molecules >= 0.5
    if has_significant_presence:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "molecules_with_Phenol",
                "operator": ">=",
                "value": 0.5,
                "scope": "proportion_of_total_molecules"
            }
        })

    has_significant_preservation = (
        total_reactions == 0 or preserved_in_reactions / total_reactions >= 0.5
    )
    if has_significant_preservation:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "preservation_of_Phenol",
                "operator": ">=",
                "value": 0.5,
                "scope": "proportion_of_total_reactions"
            }
        })

    final_product_has_directing_group = any(
        depth == 0 for _, depth in molecules_with_phenol
    )
    if final_product_has_directing_group:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Phenol",
                "position": "last_stage"
            }
        })

    if (
        has_significant_presence
        and has_significant_preservation
        and final_product_has_directing_group
    ):
        result = True
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "significant_presence_of_Phenol",
                    "significant_preservation_of_Phenol",
                    "final_product_has_Phenol"
                ]
            }
        })

    return result, findings_json
