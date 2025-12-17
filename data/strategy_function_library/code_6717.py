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


ACTIVATED_ACID_DERIVATIVES = [
    "Acyl halide",  # Includes acyl chloride, bromide, iodide
    "Anhydride",
    "Ester",
    "Primary amide",
    "Secondary amide",
    "Tertiary amide",
    "Azide",  # For acyl azide
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the activation of a carboxylic acid by its conversion into a more reactive derivative. This is identified by checking for the formation of specific functional groups defined in the ACTIVATED_ACID_DERIVATIVES list, including Acyl halide, Anhydride, Ester, Amide, and Azide.
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

    acid_activation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal acid_activation_found, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid in reactants
            acid_in_reactants = False
            for reactant in reactants:
                if checker.check_fg("Carboxylic acid", reactant):
                    acid_in_reactants = True
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                    break

            # If carboxylic acid found in reactants, check for activated derivatives in product
            if acid_in_reactants:
                # Check for activated derivatives in product that were not present in reactants
                for derivative in ACTIVATED_ACID_DERIVATIVES:
                    if checker.check_fg(derivative, product) and not checker.check_fg(
                        derivative, rsmi.split(">")[0]
                    ):
                        acid_activation_found = True
                        if derivative not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(derivative)
                        # Add the structural constraint when the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "single_reaction",
                                "description": "A reaction step must meet three conditions: 1) A 'Carboxylic acid' functional group is present in the reactants. 2) An activated acid derivative (e.g., 'Acyl halide', 'Ester') is present in the product. 3) The same activated acid derivative is NOT present in the reactants."
                            }
                        })
                        return

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return acid_activation_found, findings_json
