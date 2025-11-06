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


N_ARYLATION_AND_SNAR_REACTIONS = [
    "Nucleophilic substitution aromatic",
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "N-arylation",
    "Buchwald-Hartwig",
    "Goldberg coupling",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a late-stage nucleophilic aromatic substitution with an amine.
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

    found_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar, findings_json

        if node["type"] == "reaction" and depth <= 3:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Rely on high-fidelity reaction checkers. The original code's
                # subsequent functional group analysis was flawed and removed.
                for name in N_ARYLATION_AND_SNAR_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        found_snar = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if found_snar:
        # Add the structural constraint if any of the target reactions were found within the depth limit
        # This corresponds to the 'positional' constraint in the strategy JSON.
        # The 'targets' list in the constraint matches N_ARYLATION_AND_SNAR_REACTIONS.
        # The 'position' is 'last_4_stages' which is implicitly checked by depth <= 3.
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "targets": [
                    "Nucleophilic substitution aromatic",
                    "heteroaromatic_nuc_sub",
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                    "N-arylation",
                    "Buchwald-Hartwig",
                    "Goldberg coupling",
                    "Ullmann condensation"
                ],
                "position": "last_4_stages"
            }
        })

    return found_snar, findings_json
