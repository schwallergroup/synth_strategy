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


ESTER_FORMATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "Transesterification",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Schotten-Baumann to ester",
    "Oxidative esterification of primary alcohols",
    "Oxidation of alcohol and aldehyde to ester",
    "Acetic anhydride and alcohol to ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving late-stage ester formation. Late-stage is defined as occurring within the final 3 steps (depth 0-2). The strategy is identified if the reaction is one of the named reactions in the `ESTER_FORMATION_REACTIONS` list, or, as a fallback, if a carboxylic acid is consumed and an ester is formed, provided no ester was present in the reactants.
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

    has_late_stage_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_esterification, findings_json

        if node["type"] == "reaction" and depth <= 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                
                is_ester_formation_rxn = False
                for rxn_name in ESTER_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_ester_formation_rxn = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                if is_ester_formation_rxn:
                    has_late_stage_esterification = True
                    if {"type": "positional", "details": {"target": "ester_formation_event", "position_depth": {"operator": "<=", "value": 2}}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ester_formation_event", "position_depth": {"operator": "<=", "value": 2}}})
                else:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]
                    
                    has_acid_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", reactant):
                            has_acid_in_reactants = True
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            break

                    has_ester_in_product = checker.check_fg("Ester", product_smiles)
                    if has_ester_in_product:
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")

                    reactants_have_ester = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Ester", reactant):
                            reactants_have_ester = True
                            # No need to add Ester to functional_groups again if already added from product check
                            break

                    if has_acid_in_reactants and has_ester_in_product and not reactants_have_ester:
                        has_late_stage_esterification = True
                        if {"type": "positional", "details": {"target": "ester_formation_event", "position_depth": {"operator": "<=", "value": 2}}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ester_formation_event", "position_depth": {"operator": "<=", "value": 2}}})
                        if {"type": "co-occurrence", "details": {"targets": ["Carboxylic acid in reactants", "Ester in product"], "scope": "single_reaction_step"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Carboxylic acid in reactants", "Ester in product"], "scope": "single_reaction_step"}})
                        if {"type": "negation", "details": {"target": "Ester in reactants", "scope": "single_reaction_step"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Ester in reactants", "scope": "single_reaction_step"}})

            except Exception:
                pass

        # Continue traversing children regardless of what we found
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return has_late_stage_esterification, findings_json