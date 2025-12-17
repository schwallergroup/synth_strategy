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
    This function detects a synthetic strategy involving hydrazide formation
    in the late stage of synthesis.
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

    hydrazide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal hydrazide_formation_detected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late stage reaction (depth <= 1)
            if depth <= 1:
                # Check for the formation of an acylhydrazine group from its precursors
                product_is_acylhydrazine = checker.check_fg("Acylhydrazine", product)
                if product_is_acylhydrazine:
                    if "Acylhydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Acylhydrazine")

                reactants_have_acylhydrazine = any(checker.check_fg("Acylhydrazine", r) for r in reactants)
                if product_is_acylhydrazine and not reactants_have_acylhydrazine:
                    # Add negation constraint if met
                    if {"type": "negation", "details": {"target": "Acylhydrazine", "context": "reaction_reactants"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Acylhydrazine", "context": "reaction_reactants"}})

                    has_hydrazine = False
                    for r in reactants:
                        if checker.check_fg("Hydrazine", r):
                            has_hydrazine = True
                            if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")
                            break

                    has_carboxylic_derivative = False
                    carboxylic_derivatives_found = []
                    for r in reactants:
                        if checker.check_fg("Carboxylic acid", r):
                            has_carboxylic_derivative = True
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            carboxylic_derivatives_found.append("Carboxylic acid")
                        if checker.check_fg("Ester", r):
                            has_carboxylic_derivative = True
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                            carboxylic_derivatives_found.append("Ester")
                        if checker.check_fg("Acyl halide", r):
                            has_carboxylic_derivative = True
                            if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                            carboxylic_derivatives_found.append("Acyl halide")

                    if has_hydrazine and has_carboxylic_derivative:
                        hydrazide_formation_detected = True
                        # Add co-occurrence constraint if met
                        if {"type": "co-occurrence", "details": {"targets": ["Hydrazine", ["Carboxylic acid", "Ester", "Acyl halide"]], "context": "reaction_reactants"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Hydrazine", ["Carboxylic acid", "Ester", "Acyl halide"]], "context": "reaction_reactants"}})
                        
                        # Add positional constraint if met (late stage)
                        if {"type": "positional", "details": {"target": "hydrazide_formation", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "hydrazide_formation", "position": "late_stage"}})
                        
                        # Add named reaction if met
                        if "hydrazide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("hydrazide_formation")

        # Continue traversing with increased depth
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return hydrazide_formation_detected, findings_json
