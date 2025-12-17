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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where Boc protection is maintained
    throughout the synthesis without deprotection/reprotection steps.
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

    # Track molecules with Boc at each depth
    boc_present_at_depths = {}
    # Track if Boc protection/deprotection reactions occur
    boc_protection_reactions = set()
    boc_deprotection_reactions = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            # Check for Boc group in molecule
            mol_smiles = node["smiles"]
            has_boc = checker.check_fg("Boc", mol_smiles)

            if has_boc:
                if depth not in boc_present_at_depths:
                    boc_present_at_depths[depth] = []
                boc_present_at_depths[depth].append(mol_smiles)
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

        elif node["type"] == "reaction":
            # Extract reactants and product from reaction
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for Boc protection/deprotection reactions
            is_boc_protection = any(checker.check_reaction(r, rsmi) for r in BOC_PROTECTION_REACTIONS)
            is_boc_deprotection = any(checker.check_reaction(r, rsmi) for r in BOC_DEPROTECTION_REACTIONS)

            if is_boc_protection:
                boc_protection_reactions.add(depth)
                for r_name in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

            if is_boc_deprotection:
                boc_deprotection_reactions.add(depth)
                for r_name in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] == "mol": # From chemical to reaction node, depth increases
                new_depth = depth + 1
            # If node["type"] == "reaction", depth remains the same (reaction to chemical)
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if:
    # 1. Boc is present at multiple depths (at least 2)
    # 2. No Boc deprotection reactions occur
    # 3. At most one Boc protection reaction occurs (at the beginning)

    depths_with_boc = sorted(boc_present_at_depths.keys())
    has_multiple_depths_with_boc = len(depths_with_boc) >= 2
    if has_multiple_depths_with_boc:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "stages_with_Boc_group", "operator": ">=", "value": 2, "description": "The Boc functional group must be present in molecules at two or more different stages of the synthesis."}})

    # Check if Boc is maintained throughout synthesis
    continuous_boc = False
    if depths_with_boc:
        # Check if Boc is present in a continuous range from the target molecule
        expected_depths = list(range(min(depths_with_boc), max(depths_with_boc) + 1))
        continuous_boc = all(depth in boc_present_at_depths for depth in expected_depths)
    if continuous_boc:
        findings_json["structural_constraints"].append({"type": "sequence", "details": {"description": "The 'Boc' functional group must be present in at least one molecule at every synthesis stage between its first and last appearance in the route."}})

    # No deprotection should occur
    no_deprotection = len(boc_deprotection_reactions) == 0
    if no_deprotection:
        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "any_Boc_deprotection_reaction", "description": "No Boc deprotection reactions should occur anywhere in the synthetic route."}})

    # At most one protection reaction (typically at the beginning of synthesis)
    valid_protection = len(boc_protection_reactions) <= 1
    if valid_protection:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "Boc_protection_reaction", "operator": "<=", "value": 1, "description": "At most one Boc protection reaction can occur in the synthetic route."}})

    # Strategy is valid if all conditions are met
    strategy_present = (
        has_multiple_depths_with_boc and continuous_boc and no_deprotection and valid_protection
    )

    return strategy_present, findings_json
