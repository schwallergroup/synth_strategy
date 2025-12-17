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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Ester with ammonia to amide",
    "Carboxylic acid to amide conversion",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis route uses a late-stage amide formation strategy,
    where an amide bond is formed in the final step of the synthesis.
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

    final_step_is_amide_formation = False
    min_depth_reactions = {}  # Track minimum depth for each reaction

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_formation, findings_json

        if node["type"] == "reaction":
            rxn_id = node["metadata"].get("reaction_hash", node["metadata"].get("rsmi", "unknown"))
            if rxn_id not in min_depth_reactions or depth < min_depth_reactions[rxn_id]:
                min_depth_reactions[rxn_id] = depth

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

        if node["type"] == "reaction":
            if depth <= 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                total_amide_count_product = sum([
                    len(checker.get_fg_atom_indices("Primary amide", product_smiles)),
                    len(checker.get_fg_atom_indices("Secondary amide", product_smiles)),
                    len(checker.get_fg_atom_indices("Tertiary amide", product_smiles)),
                ])
                if len(checker.get_fg_atom_indices("Primary amide", product_smiles)) > 0 and "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if len(checker.get_fg_atom_indices("Secondary amide", product_smiles)) > 0 and "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if len(checker.get_fg_atom_indices("Tertiary amide", product_smiles)) > 0 and "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                total_amide_count_reactants = sum(
                    len(checker.get_fg_atom_indices("Primary amide", r)) +
                    len(checker.get_fg_atom_indices("Secondary amide", r)) +
                    len(checker.get_fg_atom_indices("Tertiary amide", r))
                    for r in reactants_smiles
                )

                is_amide_formation_reaction = False
                for name in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_amide_formation_reaction = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                        break

                if is_amide_formation_reaction and (total_amide_count_product > total_amide_count_reactants):
                    final_step_is_amide_formation = True
                    # Add structural constraint if the condition is met
                    if {"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}})

    dfs_traverse(route)

    return final_step_is_amide_formation, findings_json
