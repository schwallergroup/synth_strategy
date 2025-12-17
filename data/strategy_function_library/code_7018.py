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


TRIAZOLE_FORMING_REACTIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen Cu-catalyzed 1,4-subst",
    "Huisgen Ru-catalyzed 1,5_subst",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Huisgen disubst-alkyne",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving a triazole-containing fragment, identified by either the presence of a triazole ring in any molecule or the use of a specific triazole-forming reaction. This function specifically checks for Huisgen cycloaddition reactions, including: 'Huisgen alkyne-azide 1,3 dipolar cycloaddition', 'Huisgen Cu-catalyzed 1,4-subst', 'Huisgen Ru-catalyzed 1,5_subst', 'Huisgen 1,3 dipolar cycloaddition', 'Huisgen alkene-azide 1,3 dipolar cycloaddition', and 'Huisgen disubst-alkyne'.
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

    has_triazole_fragment = False
    has_triazole_forming_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_triazole_fragment, has_triazole_forming_reaction, findings_json

        if node["type"] == "mol":
            smiles = node["smiles"]
            if checker.check_ring("triazole", smiles):
                has_triazole_fragment = True
                if "triazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("triazole")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            for name in TRIAZOLE_FORMING_REACTIONS:
                if checker.check_reaction(name, rxn_smiles):
                    has_triazole_forming_reaction = True
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    # No break here, as multiple reactions might match, and we want to record all of them

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same

            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = has_triazole_fragment or has_triazole_forming_reaction
    return result, findings_json
