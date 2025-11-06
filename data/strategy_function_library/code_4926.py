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


BORYLATION_REACTION_TYPES = [
    "Preparation of boronic acids",
    "Preparation of boronic acids without boronic ether",
    "Preparation of boronic ethers",
    "Preparation of boronic acids from trifluoroborates",
    "Synthesis of boronic acids",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves a late-stage borylation (C-B bond formation)
    in the final steps of the synthesis (final 3 steps).
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

    found_borylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_borylation, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Late stage: final 3 steps
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if the reaction is a known borylation reaction
                for rxn_type in BORYLATION_REACTION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        found_borylation = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        if {"type": "positional", "details": {"target": "borylation", "position": "late_stage", "description": "The borylation event must occur within the final 3 reaction steps (depth <= 2)."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "borylation", "position": "late_stage", "description": "The borylation event must occur within the final 3 reaction steps (depth <= 2)."}})
                        return

                # Check if product contains boron-containing functional groups
                product_has_boron = False
                if checker.check_fg("Boronic acid", product_smiles):
                    product_has_boron = True
                    if "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                if checker.check_fg("Boronic ester", product_smiles):
                    product_has_boron = True
                    if "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")

                if not product_has_boron:
                    return

                # Check if reactants lack boron-containing functional groups
                reactants_have_boron = False
                for reactant in reactants_smiles:
                    if checker.check_fg("Boronic acid", reactant):
                        reactants_have_boron = True
                        if "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                    if checker.check_fg("Boronic ester", reactant):
                        reactants_have_boron = True
                        if "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")

                # True borylation: product has boron but reactants don't
                if product_has_boron and not reactants_have_boron:
                    found_borylation = True
                    if "borylation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("borylation")
                    if {"type": "positional", "details": {"target": "borylation", "position": "late_stage", "description": "The borylation event must occur within the final 3 reaction steps (depth <= 2)."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "borylation", "position": "late_stage", "description": "The borylation event must occur within the final 3 reaction steps (depth <= 2)."}})
            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return found_borylation, findings_json
