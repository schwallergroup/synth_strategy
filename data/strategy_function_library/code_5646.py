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


SNAR_REACTIONS_OF_INTEREST = [
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "heteroaromatic_nuc_sub",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects specific nucleophilic aromatic substitution (SNAr) reactions by checking for the presence of the following named reactions: nucl_sub_aromatic_ortho_nitro, nucl_sub_aromatic_para_nitro, and heteroaromatic_nuc_sub.
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

    has_nucleophilic_substitution = False

    def is_nucleophilic_substitution(reaction_smiles):
        # Check for specific SNAr reactions from the defined list
        nonlocal findings_json
        found_any_snar = False
        for reaction_name in SNAR_REACTIONS_OF_INTEREST:
            if checker.check_reaction(reaction_name, reaction_smiles):
                if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                found_any_snar = True
        return found_any_snar

    def dfs_traverse(node, depth):
        nonlocal has_nucleophilic_substitution

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                reaction_smiles = node["metadata"]["mapped_reaction_smiles"]

                # Check for nucleophilic aromatic substitution
                if is_nucleophilic_substitution(reaction_smiles):
                    has_nucleophilic_substitution = True

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route, 0) # Initial depth is 0

    return has_nucleophilic_substitution, findings_json
