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
    This function detects a synthetic strategy involving late-stage amide formation
    preceded by ester hydrolysis and early ether formation.
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

    # Initialize flags to track strategy components
    has_amide_formation = False
    has_ester_hydrolysis = False
    has_ether_formation = False

    # Track the depths at which each reaction occurs (initialize with large values)
    amide_depth = float("inf")
    hydrolysis_depth = float("inf")
    ether_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal has_amide_formation, has_ester_hydrolysis, has_ether_formation
        nonlocal amide_depth, hydrolysis_depth, ether_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for amide formation reactions
            amide_formation_reactions = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Schotten-Baumann_amide",
                "Carboxylic acid to amide conversion",
                "Acylation of primary amines",
                "Acylation of secondary amines",
            ]

            for reaction_name in amide_formation_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    has_amide_formation = True
                    amide_depth = min(amide_depth, depth)
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            # Check for ester hydrolysis
            ester_hydrolysis_reactions = [
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Ester saponification (methyl deprotection)",
                "Ester saponification (alkyl deprotection)",
                "COOH ethyl deprotection",
                "Carboxylate to carboxylic acid",
            ]

            for reaction_name in ester_hydrolysis_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    has_ester_hydrolysis = True
                    hydrolysis_depth = min(hydrolysis_depth, depth)
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            # Check for ether formation
            ether_formation_reactions = [
                "Williamson Ether Synthesis",
                "{Williamson ether}",
                "Alcohol to ether",
                "Mitsunobu aryl ether",
                "Williamson Ether Synthesis (intra to epoxy)",
                "Chan-Lam etherification",
                "Ullmann-Goldberg Substitution aryl alcohol",
            ]

            for reaction_name in ether_formation_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    has_ether_formation = True
                    ether_depth = min(ether_depth, depth)
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the complete strategy is present with correct sequence
    late_stage_amide = amide_depth <= 1  # Amide formation is late-stage (depth 0 or 1)

    # In retrosynthesis, we traverse from product to reactants, so the sequence should be:
    # amide formation (lowest depth) -> ester hydrolysis -> ether formation (highest depth)
    correct_sequence = (
        (amide_depth < hydrolysis_depth < ether_depth)
        if (has_amide_formation and has_ester_hydrolysis and has_ether_formation)
        else False
    )

    strategy_present = (
        has_amide_formation
        and has_ester_hydrolysis
        and has_ether_formation
        and late_stage_amide
        and correct_sequence
    )

    # Record structural constraints if met
    if has_amide_formation and has_ester_hydrolysis and has_ether_formation:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "amide_formation",
                    "ester_hydrolysis",
                    "ether_formation"
                ]
            }
        })
    
    if has_amide_formation and late_stage_amide:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_formation",
                "position": "depth <= 1"
            }
        })

    if correct_sequence:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "ether_formation",
                    "ester_hydrolysis",
                    "amide_formation"
                ],
                "description": "The order reflects the forward synthesis direction, inferred from the retrosynthesis depth check where depth(amide) < depth(hydrolysis) < depth(ether)."
            }
        })

    return strategy_present, findings_json
