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
    Detects a strategy involving protection/deprotection steps combined with
    ring modification (opening, reduction, or aromatization).
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

    # Initialize tracking variables
    has_protection = False
    has_deprotection = False
    has_ring_modification = False

    def dfs_traverse(node, depth):
        nonlocal has_protection, has_deprotection, has_ring_modification, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for protection reactions
                protection_reactions = [
                    "Alcohol protection with silyl ethers",
                    "Boc amine protection",
                    "Protection of carboxylic acid",
                    "Aldehyde or ketone acetalization",
                    "Esterification of Carboxylic Acids",
                    "Schotten-Baumann to ester",
                    "Formation of Sulfonic Esters"
                ]
                for r_name in protection_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        has_protection = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for deprotection reactions
                deprotection_reactions = [
                    "Alcohol deprotection from silyl ethers",
                    "Boc amine deprotection",
                    "Deprotection of carboxylic acid",
                    "Acetal hydrolysis to aldehyde",
                    "Ketal hydrolysis to ketone",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)"
                ]
                for r_name in deprotection_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        has_deprotection = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for specific ring-modifying reactions
                ring_modification_reactions = [
                    "Diels-Alder",
                    "Retro-Diels-Alder from oxazole",
                    "Arene hydrogenation"
                ]
                for r_name in ring_modification_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        has_ring_modification = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases for reaction children
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route, 0) # Initial call with depth 0

    # Strategy is present if we have both protection/deprotection and ring modification
    strategy_present = (has_protection or has_deprotection) and has_ring_modification

    if strategy_present:
        # Add the structural constraint if the strategy is detected
        findings_json["structural_constraints"].append(
            {
                "type": "co-occurrence",
                "details": {
                    "description": "Requires at least one reaction from the protection/deprotection group to co-occur with at least one reaction from the ring modification group anywhere in the synthesis route.",
                    "targets": [
                        [
                            "Alcohol protection with silyl ethers",
                            "Boc amine protection",
                            "Protection of carboxylic acid",
                            "Aldehyde or ketone acetalization",
                            "Esterification of Carboxylic Acids",
                            "Schotten-Baumann to ester",
                            "Formation of Sulfonic Esters",
                            "Alcohol deprotection from silyl ethers",
                            "Boc amine deprotection",
                            "Deprotection of carboxylic acid",
                            "Acetal hydrolysis to aldehyde",
                            "Ketal hydrolysis to ketone",
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                            "Ester saponification (methyl deprotection)",
                            "Ester saponification (alkyl deprotection)"
                        ],
                        [
                            "Diels-Alder",
                            "Retro-Diels-Alder from oxazole",
                            "Arene hydrogenation"
                        ]
                    ]
                }
            }
        )

    return strategy_present, findings_json
