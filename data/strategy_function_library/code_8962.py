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
    This function detects late-stage ketone reduction (C=O to C-OH).
    Late stage is defined as occurring at depth 0 or 1.

    In retrosynthetic analysis, we look for the reverse transformation:
    secondary alcohol in reactants being converted to ketone in product.
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

    late_stage_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_reduction_found, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for alcohol oxidation in retrosynthesis
                oxidation_reaction_1_name = "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones"
                oxidation_reaction = checker.check_reaction(
                    oxidation_reaction_1_name, rsmi
                )
                if oxidation_reaction:
                    if oxidation_reaction_1_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(oxidation_reaction_1_name)

                if not oxidation_reaction:
                    oxidation_reaction_2_name = "Oxidation of secondary alcohol to ketone"
                    oxidation_reaction = checker.check_reaction(
                        oxidation_reaction_2_name, rsmi
                    )
                    if oxidation_reaction:
                        if oxidation_reaction_2_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(oxidation_reaction_2_name)

                if oxidation_reaction:
                    # In retrosynthesis: alcohol in reactants, ketone in product
                    has_secondary_alcohol = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Secondary alcohol", reactant):
                            has_secondary_alcohol = True
                            if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                            break

                    has_ketone = checker.check_fg("Ketone", product_smiles)
                    if has_ketone:
                        if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")

                    if has_ketone and has_secondary_alcohol:
                        late_stage_reduction_found = True
                        # Record structural constraint: co-occurrence
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "alcohol_oxidation",
                                    "Secondary alcohol",
                                    "Ketone"
                                ],
                                "scope": "single_reaction_step"
                            }
                        })
                        # Record structural constraint: positional (late_stage)
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "alcohol_oxidation",
                                "position": "late_stage",
                                "constraint": "depth <= 1"
                            }
                        })

            except Exception:
                pass

        # Process children
        for child in node.get("children", []):
            if late_stage_reduction_found:
                break
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return late_stage_reduction_found, findings_json
