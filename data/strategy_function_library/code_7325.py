from typing import Tuple, Dict, List
import copy
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
    Detects if a sulfonamide group is preserved across a synthetic route. This is true if the group is present in at least three non-starting material molecules and is never lost in any reaction step where it was present in a major reactant.
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

    # Track if sulfonamide is preserved through each reaction
    preserved_through_all = True

    # Track molecules with sulfonamide for reporting
    molecules_with_sulfonamide = 0

    def dfs_traverse(node, depth=0):
        nonlocal preserved_through_all, molecules_with_sulfonamide, findings_json

        # Process molecule nodes (excluding starting materials)
        if node["type"] == "mol" and not node.get("in_stock", False):
            # Skip very small molecules (like water, HBr, etc.)
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.GetNumAtoms() > 2:  # Skip small molecules
                has_sulfonamide = checker.check_fg("Sulfonamide", node["smiles"])

                if has_sulfonamide:
                    molecules_with_sulfonamide += 1
                    if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

        # Process reaction nodes to check if sulfonamide is preserved
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if sulfonamide exists in product and reactants
            product_has_sulfonamide = checker.check_fg("Sulfonamide", product)

            # Filter out small molecules from reactants
            significant_reactants = []
            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol and mol.GetNumAtoms() > 2:
                    significant_reactants.append(r)

            # Check if any significant reactant has sulfonamide
            reactant_has_sulfonamide = any(
                checker.check_fg("Sulfonamide", r) for r in significant_reactants
            )

            # If a reactant has sulfonamide but product doesn't, sulfonamide was lost
            if reactant_has_sulfonamide and not product_has_sulfonamide:
                preserved_through_all = False

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Final result calculation
    result = molecules_with_sulfonamide >= 3 and preserved_through_all

    # Add structural constraints to findings_json based on the final result
    if molecules_with_sulfonamide >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "molecule_with_Sulfonamide",
                "operator": ">=",
                "value": 3
            }
        })
    
    if preserved_through_all:
        # This constraint is met if sulfonamide was NOT lost in any reaction
        # The original JSON describes a 'negation' constraint, meaning the absence of 'reaction_losing_Sulfonamide'
        # If preserved_through_all is True, then 'reaction_losing_Sulfonamide' did not occur.
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "reaction_losing_Sulfonamide"
            }
        })

    # Return True if at least 3 non-starting molecules have a sulfonamide group
    # AND sulfonamide is preserved through all reactions
    return result, findings_json
