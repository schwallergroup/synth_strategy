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


BORYLATION_REACTIONS = [
    "Preparation of boronic acids",
    "Preparation of boronic acids without boronic ether",
    "Preparation of boronic ethers",
    "Preparation of boronic acids from trifluoroborates",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of boronic acids or esters by identifying reactions from a defined list of borylation methods, such as 'Preparation of boronic acids' and 'Preparation of boronic ethers'.
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

    has_borylation = False

    def dfs_traverse(node):
        nonlocal has_borylation, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a borylation reaction using the reaction checker
            for r_name in BORYLATION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    print(f"Found borylation reaction by reaction type: {rsmi}")
                    has_borylation = True
                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    # If multiple borylation reactions can be found, remove 'break' to record all.
                    # For this specific problem, we can break after the first match for efficiency
                    # if only existence matters for the boolean flag.
                    # However, to record all, we should not break.
                    # For now, keeping the original logic's implication of 'any' match.
                    break # Break after finding the first matching reaction to set has_borylation and record it.

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_borylation, findings_json
