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
    This function detects if nitro groups are preserved throughout the synthesis.
    It checks if nitro groups are present in the final product and maintained
    through the synthesis route.
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

    # Track molecules with nitro groups and their depth
    mols_with_nitro = []
    max_depth = 0
    nitro_reactions_preserved = 0
    total_nitro_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal mols_with_nitro, max_depth, nitro_reactions_preserved, total_nitro_reactions, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Skip starting materials
            if not node.get("in_stock", False):
                # Check for nitro group using the checker function
                if checker.check_fg("Nitro group", mol_smiles):
                    mols_with_nitro.append((mol_smiles, depth))
                    if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

        elif node["type"] == "reaction":
            # Check if nitro groups are preserved through this reaction
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has nitro group
                product_has_nitro = checker.check_fg("Nitro group", product)

                # Check if any reactant has nitro group
                reactant_has_nitro = any(checker.check_fg("Nitro group", r) for r in reactants)

                if reactant_has_nitro:
                    total_nitro_reactions += 1
                    if product_has_nitro:
                        nitro_reactions_preserved += 1

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check conditions for nitro group preservation
    result = False

    if not mols_with_nitro:
        return result, findings_json

    # Check if final product (depth 0) has nitro group
    final_product_has_nitro = any(depth == 0 for _, depth in mols_with_nitro)
    if final_product_has_nitro:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitro group", "position": "last_stage"}})

    # Check if nitro groups are present in early stages
    # Note: depth=max_depth is the first step, depth=0 is the last.
    # This checks if a nitro group appears in the first 30% of steps.
    early_stage_has_nitro = any(depth >= max_depth * 0.7 for _, depth in mols_with_nitro)
    if early_stage_has_nitro:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitro group", "position": "first_30_percent_of_route"}})

    # Calculate percentage of reactions that preserved nitro groups
    nitro_reaction_preservation = (
        nitro_reactions_preserved / total_nitro_reactions if total_nitro_reactions > 0 else 1.0
    )

    # Return True if:
    # 1. Final product has a nitro group.
    # 2. A nitro group is present in an early stage of the synthesis.
    # 3. At least 80% of reactions involving a nitro group preserve it.
    result = (
        final_product_has_nitro
        and early_stage_has_nitro
        and nitro_reaction_preservation >= 0.8
    )

    if result:
        # If the overall condition is met, it implies at least one nitro group was found.
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "Nitro group", "operator": ">=", "value": 1}})
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "ratio_of_nitro_preserving_reactions", "operator": ">=", "value": 0.8}})

    return result, findings_json