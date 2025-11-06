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


CROSS_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
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
    "Heck terminal vinyl",
    "Heck_non-terminal_vinyl",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
    "Kumada cross-coupling",
    "Hiyama-Denmark Coupling",
    "Aryllithium cross-coupling",
    "Buchwald-Hartwig",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the route contains a cross-coupling reaction by checking against a predefined list of named reactions.
    This list includes variants of Suzuki, Stille, Heck, Sonogashira, Negishi, Kumada, and Buchwald-Hartwig reactions.
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

    found_transformation = False

    def dfs_traverse(node):
        nonlocal found_transformation, findings_json

        if found_transformation:
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for various cross-coupling reactions directly
            for reaction_type in CROSS_COUPLING_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    found_transformation = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    return

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_transformation, findings_json
