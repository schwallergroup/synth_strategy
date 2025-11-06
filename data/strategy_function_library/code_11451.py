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


LATE_STAGE_ORGANOMETALLIC_REACTIONS = [
    # Grignard reactions
    "Grignard from aldehyde to alcohol",
    "Grignard from ketone to alcohol",
    "Grignard from nitrile to ketone",
    "Grignard with CO2 to carboxylic acid",
    "Olefination of ketones with Grignard reagents",
    "Olefination of aldehydes with Grignard reagents",
    "Preparation of trialkylsilanes with Grignard reagents",
    "Formation of Grignard reagents",
    "Grignard_carbonyl",
    "Grignard_alcohol",
    # Cross-coupling reactions using organometallics
    "Kumada cross-coupling",
    "Negishi coupling",
    "Stille reaction_vinyl",
    "Stille reaction_aryl",
    "Stille reaction_benzyl",
    "Stille reaction_allyl",
    "Stille reaction_vinyl OTf",
    "Stille reaction_aryl OTf",
    "Stille reaction_benzyl OTf",
    "Stille reaction_allyl OTf",
    "Stille reaction_other",
    "Stille reaction_other OTf",
    "Hiyama-Denmark Coupling",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis employs a late-stage (final or penultimate step) functionalization using a specific set of organometallic-mediated reactions. The check includes various Grignard additions and a comprehensive list of cross-coupling reactions such as Kumada, Negishi, Stille, and Hiyama-Denmark.
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

    found_late_organometallic = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_organometallic, findings_json

        if node["type"] == "reaction":
            # Extract reaction SMILES
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a late-stage reaction (depth 0 or 1)
                node_depth = node["metadata"].get("depth", depth)
                if node_depth <= 1:
                    # Check for specific organometallic reaction types
                    for reaction_type in LATE_STAGE_ORGANOMETALLIC_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            found_late_organometallic = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            # No break here, as we want to record all found reactions

        # Continue traversing with updated depth
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # Current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if found_late_organometallic:
        # Add the structural constraint if any late-stage organometallic reaction was found
        findings_json["structural_constraints"].append(
            {
                "type": "positional",
                "details": {
                    "targets": [
                        "Grignard from aldehyde to alcohol",
                        "Grignard from ketone to alcohol",
                        "Grignard from nitrile to ketone",
                        "Grignard with CO2 to carboxylic acid",
                        "Olefination of ketones with Grignard reagents",
                        "Olefination of aldehydes with Grignard reagents",
                        "Preparation of trialkylsilanes with Grignard reagents",
                        "Formation of Grignard reagents",
                        "Grignard_carbonyl",
                        "Grignard_alcohol",
                        "Kumada cross-coupling",
                        "Negishi coupling",
                        "Stille reaction_vinyl",
                        "Stille reaction_aryl",
                        "Stille reaction_benzyl",
                        "Stille reaction_allyl",
                        "Stille reaction_vinyl OTf",
                        "Stille reaction_aryl OTf",
                        "Stille reaction_benzyl OTf",
                        "Stille reaction_allyl OTf",
                        "Stille reaction_other",
                        "Stille reaction_other OTf",
                        "Hiyama-Denmark Coupling",
                        "Aryllithium cross-coupling"
                    ],
                    "position": [
                        "last_stage",
                        "penultimate_stage"
                    ]
                }
            }
        )

    return found_late_organometallic, findings_json
