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


SULFONAMIDE_FORMATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage (final two steps) installation of or coupling with a sulfonamide. The strategy is flagged if a sulfonamide is formed via a known reaction (as defined in the `SULFONAMIDE_FORMATION_REACTIONS` list) or if a reaction couples a reactant that already contains a sulfonamide moiety.
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

    found_coupling = False

    # Define the structural constraints from the strategy JSON for easy lookup
    # In a real scenario, these would be passed in or loaded dynamically.
    # For this exercise, we hardcode them based on the provided JSON.
    sulfonamide_installation_constraint = {
      "type": "co-occurrence",
      "details": {
        "description": "A named sulfonamide synthesis reaction occurs within the last 3 stages of the synthesis, and the product of that reaction contains a sulfonamide group. This represents the 'installation' pathway.",
        "scope": "single_reaction_step",
        "targets": [
          "position_within_last_3_stages",
          "named_sulfonamide_synthesis",
          "product_contains_sulfonamide"
        ]
      }
    }
    sulfonamide_coupling_constraint = {
      "type": "co-occurrence",
      "details": {
        "description": "A reaction within the last 3 stages involves a reactant that contains a sulfonamide, and the product of that reaction also contains a sulfonamide. This represents the 'coupling' pathway.",
        "scope": "single_reaction_step",
        "targets": [
          "position_within_last_3_stages",
          "reactant_contains_sulfonamide",
          "product_contains_sulfonamide"
        ]
      }
    }

    def dfs_traverse(node, depth=0):
        nonlocal found_coupling, findings_json

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Late stage (up to 2 steps from final product)
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                sulfonamide_in_product = checker.check_fg("Sulfonamide", product)
                if sulfonamide_in_product and "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

                is_sulfonamide_rxn = False
                for rxn_name in SULFONAMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_sulfonamide_rxn = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                sulfonamide_in_reactants = False
                for reactant in reactants:
                    if checker.check_fg("Sulfonamide", reactant):
                        sulfonamide_in_reactants = True
                        if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                        break

                if is_sulfonamide_rxn and sulfonamide_in_product:
                    found_coupling = True
                    # Add the installation constraint if not already present
                    if sulfonamide_installation_constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(sulfonamide_installation_constraint)

                if sulfonamide_in_reactants and sulfonamide_in_product:
                    found_coupling = True
                    # Add the coupling constraint if not already present
                    if sulfonamide_coupling_constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(sulfonamide_coupling_constraint)

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_coupling, findings_json