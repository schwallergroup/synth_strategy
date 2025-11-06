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


PHENOL_ETHER_FORMATION_REACTIONS = [
    "Williamson Ether Synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Mitsunobu aryl ether",
    "Chan-Lam etherification",
    "Ullmann-Goldberg Substitution aryl alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the formation of a phenol ether in the early stages of synthesis (depth >= 2) via a specific set of named reactions. This strategy identifies reactions where a phenol is a reactant, an aryl ether is a product, and the transformation is confirmed to be one of the following types: Williamson Ether Synthesis, Williamson Ether Synthesis (intra to epoxy), Mitsunobu aryl ether, Chan-Lam etherification, or Ullmann-Goldberg Substitution aryl alcohol.
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

    phenol_disconnection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal phenol_disconnection_found, findings_json

        if node["type"] == "reaction" and depth >= 2:  # Early stage (higher depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for phenol in reactants
                phenol_reactant = None
                for r_smiles in reactants_smiles:
                    if checker.check_fg("Phenol", r_smiles):
                        phenol_reactant = r_smiles
                        if "Phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                        break

                if phenol_reactant:
                    # Check for aryl ether in product (more specific and robust)
                    if checker.check_fg("Aryl ether", product_smiles):
                        if "Aryl ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aryl ether")

                        # Check if this is a known ether formation reaction
                        reaction_found = False
                        for rxn in PHENOL_ETHER_FORMATION_REACTIONS:
                            if checker.check_reaction(rxn, rsmi):
                                if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                                reaction_found = True
                                break
                        
                        if reaction_found:
                            phenol_disconnection_found = True
                            # Add structural constraints if all conditions are met
                            if {"type": "positional", "details": {"target": "phenol_ether_formation", "position": "early_stage"}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "phenol_ether_formation", "position": "early_stage"}})
                            if {"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["Phenol", "Aryl ether", "any_phenol_ether_formation_reaction"]}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["Phenol", "Aryl ether", "any_phenol_ether_formation_reaction"]}})

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return phenol_disconnection_found, findings_json
