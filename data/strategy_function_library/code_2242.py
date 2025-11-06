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


THIOUREA_DERIVED_HETEROCYCLE_REACTIONS = [
    "benzothiazole",
    "benzimidazole_derivatives_aldehyde",
    "benzoxazole",
    "thiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the use of a thiourea-containing reactant as a precursor for the formation of specific heterocycles, including benzothiazole, benzimidazole, benzoxazole, and thiazole.
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

    # Track if we found the strategy
    found_thiourea_intermediate = False
    thiourea_used_for_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal found_thiourea_intermediate, thiourea_used_for_heterocycle, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant has thiourea
            has_thiourea_in_reactants = False
            for reactant in reactants_smiles:
                if checker.check_fg("Thiourea", reactant):
                    has_thiourea_in_reactants = True
                    found_thiourea_intermediate = True
                    if "Thiourea" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Thiourea")
                    break

            # If thiourea found in reactants, check if it's used for heterocycle formation
            if has_thiourea_in_reactants:
                # Check if thiourea is not present in the product (consumed)
                thiourea_consumed = not checker.check_fg("Thiourea", product_smiles)
                if thiourea_consumed:
                    # This corresponds to the negation structural constraint
                    if {"type": "negation", "details": {"target": "Thiourea", "comment": "The Thiourea functional group must be consumed in the key reaction step and not be present in the product of that step."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Thiourea", "comment": "The Thiourea functional group must be consumed in the key reaction step and not be present in the product of that step."}})

                # Check if this is a heterocycle formation reaction
                is_heterocycle_formation = False
                for rxn_type in THIOUREA_DERIVED_HETEROCYCLE_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_heterocycle_formation = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        # Also add "ring_formation" if it's a heterocycle reaction
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        break

                # If thiourea is consumed and it's a known heterocycle formation reaction, it's a match.
                if thiourea_consumed and is_heterocycle_formation:
                    thiourea_used_for_heterocycle = True
                    # This corresponds to the co-occurrence structural constraint
                    if {"type": "co-occurrence", "details": {"targets": ["Thiourea", "ring_formation"], "comment": "A reaction step must involve a Thiourea functional group and result in a ring formation (e.g., benzothiazole, thiazole)."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Thiourea", "ring_formation"], "comment": "A reaction step must involve a Thiourea functional group and result in a ring formation (e.g., benzothiazole, thiazole)."}})

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if thiourea was used as a precursor for heterocycle formation
    result = found_thiourea_intermediate and thiourea_used_for_heterocycle
    return result, findings_json
