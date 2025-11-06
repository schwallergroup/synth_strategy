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
    Detects if the synthesis involves diazonium chemistry by identifying either the explicit formation or presence of a diazonium salt intermediate. It flags a route if a molecule contains a diazonium group (via SMILES string matching) or if a reaction is classified as a 'Diazo addition'.
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

    has_diazonium_chemistry = False

    # The original logic for tracking precursors was removed to prevent false positives.
    # The variable 'potential_diazonium_precursors' is no longer used.
    potential_diazonium_precursors = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_diazonium_chemistry, findings_json

        if has_diazonium_chemistry:
            return

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this molecule contains a diazonium group
            if "N#N" in mol_smiles or "N=[N+]" in mol_smiles or "[N+]#N" in mol_smiles:
                has_diazonium_chemistry = True
                if "diazonium" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("diazonium")
                return

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]

            # Check if this is a known diazo addition reaction
            if checker.check_reaction("Diazo addition", rsmi):
                has_diazonium_chemistry = True
                if "Diazo addition" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Diazo addition")
                return

            # Check for diazonium intermediates in product
            if "N#N" in product or "N=[N+]" in product or "[N+]#N" in product:
                has_diazonium_chemistry = True
                if "diazonium" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("diazonium")
                return

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol', depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return has_diazonium_chemistry, findings_json
