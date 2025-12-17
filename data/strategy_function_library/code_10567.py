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
    This function detects a synthetic strategy involving attachment of a dicyanoethylene
    group to a thiazole ring, typically through Knoevenagel condensation of an aldehyde.
    """
    strategy_detected = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # In a real implementation, max_depth would be calculated from the route object.
    # We assume it's available for the traversal.
    max_depth = 10  # Placeholder

    def dfs_traverse(node, depth, max_depth):
        nonlocal strategy_detected, findings_json

        if strategy_detected:
            return  # Early return if strategy already detected

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains a thiazole ring
                if checker.check_ring("thiazole", product_smiles):
                    findings_json["atomic_checks"]["ring_systems"].append("thiazole")
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    # Correct SMARTS for a dicyanoethylene group
                    dicyanoethylene_pattern = Chem.MolFromSmarts("[#6]=[#6]([#6]#[#7])[#6]#[#7]")

                    # Check if the product has the dicyanoethylene group
                    if product_mol and product_mol.HasSubstructMatch(dicyanoethylene_pattern):
                        findings_json["atomic_checks"]["functional_groups"].append("dicyanoethylene")
                        # Check reactants for a thiazole-aldehyde precursor
                        for reactant_smiles in reactants_smiles:
                            if checker.check_ring("thiazole", reactant_smiles):
                                if "thiazole" not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append("thiazole")
                            if checker.check_fg("Aldehyde", reactant_smiles):
                                findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                                # Confirm the dicyanoethylene group was FORMED in this step
                                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                                if reactant_mol and not reactant_mol.HasSubstructMatch(dicyanoethylene_pattern):
                                    strategy_detected = True
                                    findings_json["atomic_checks"]["named_reactions"].append("dicyanoethylene_formation")
                                    findings_json["structural_constraints"].append({
                                        "type": "co-occurrence",
                                        "details": {
                                            "targets": [
                                                "thiazole",
                                                "Aldehyde",
                                                "dicyanoethylene_formation"
                                            ],
                                            "scope": "single_reaction"
                                        }
                                    })
                                    return  # Strategy found, exit
            except Exception:
                # Silently ignore errors in SMILES processing for robustness
                pass

        # Process children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth, max_depth)

    # Start traversal from the root
    dfs_traverse(route, 1, max_depth)
    return strategy_detected, findings_json
