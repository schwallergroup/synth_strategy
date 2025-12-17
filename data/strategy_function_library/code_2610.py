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
    Detects if the synthesis route involves reduction of a nitro group followed by cyclization.
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

    nitro_reduction = False
    cyclization_step_found = False
    nitro_reduction_step_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction, cyclization_step_found, nitro_reduction_step_found, findings_json

        if (
            node["type"] == "reaction" and depth <= 3
        ):  # Expanded to depth 3 to catch more potential patterns
            try:
                # Get reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for heterocycle formation reactions
                is_heterocycle_formation = checker.check_reaction(
                    "Formation of NOS Heterocycles", rsmi
                )
                if is_heterocycle_formation:
                    findings_json["atomic_checks"]["named_reactions"].append("Formation of NOS Heterocycles")

                is_indole_formation = checker.check_reaction("Fischer indole", rsmi)
                if is_indole_formation:
                    findings_json["atomic_checks"]["named_reactions"].append("Fischer indole")

                is_benzimidazole_formation_aldehyde = checker.check_reaction(
                    "benzimidazole_derivatives_aldehyde", rsmi
                )
                if is_benzimidazole_formation_aldehyde:
                    findings_json["atomic_checks"]["named_reactions"].append("benzimidazole_derivatives_aldehyde")

                is_benzimidazole_formation_ester = checker.check_reaction(
                    "benzimidazole_derivatives_carboxylic-acid/ester", rsmi
                )
                if is_benzimidazole_formation_ester:
                    findings_json["atomic_checks"]["named_reactions"].append("benzimidazole_derivatives_carboxylic-acid/ester")

                # Case 4: Two-step process - cyclization following nitro reduction
                if depth <= 2:  # Only check for two-step process in late stages
                    # Check if this is a cyclization step
                    if (
                        is_heterocycle_formation
                        or is_indole_formation
                        or is_benzimidazole_formation_aldehyde
                        or is_benzimidazole_formation_ester
                    ):
                        cyclization_step_found = True
                        # Record positional constraint if cyclization is found at depth <= 2
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": [
                                    "Formation of NOS Heterocycles",
                                    "Fischer indole",
                                    "benzimidazole_derivatives_aldehyde",
                                    "benzimidazole_derivatives_carboxylic-acid/ester"
                                ],
                                "position": "depth <= 2"
                            }
                        })

                        # Check all children for nitro reduction
                        for child in node.get("children", []):
                            if child["type"] == "reaction":
                                try:
                                    child_rsmi = child["metadata"]["mapped_reaction_smiles"]

                                    # Check for explicit nitro reduction
                                    if checker.check_reaction(
                                        "Reduction of nitro groups to amines", child_rsmi
                                    ):
                                        nitro_reduction_step_found = True
                                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                                        # If both cyclization and nitro reduction are found, set the main flag
                                        nitro_reduction = True
                                        # Record sequence constraint
                                        findings_json["structural_constraints"].append({
                                            "type": "sequence",
                                            "details": {
                                                "event_A": [
                                                    "Formation of NOS Heterocycles",
                                                    "Fischer indole",
                                                    "benzimidazole_derivatives_aldehyde",
                                                    "benzimidazole_derivatives_carboxylic-acid/ester"
                                                ],
                                                "event_B": [
                                                    "Reduction of nitro groups to amines"
                                                ],
                                                "relationship": "B_is_direct_precursor_of_A"
                                            }
                                        })
                                        break # Found the nitro reduction, no need to check other children
                                except Exception:
                                    pass
            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return nitro_reduction, findings_json