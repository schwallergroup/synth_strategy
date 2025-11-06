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


C_N_COUPLING_REACTIONS = [
    "Ullmann-Goldberg Substitution amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Buchwald-Hartwig",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies specific C-N cross-coupling reactions, such as Buchwald-Hartwig and
    Ullmann-type couplings, by checking against a list of named reactions.
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

    # Initialize tracking variable
    has_aromatic_nucleophilic_substitution = False

    def dfs_traverse(node, depth, max_depth):
        nonlocal has_aromatic_nucleophilic_substitution, findings_json

        if node["type"] == "reaction":
            # Extract reactants and products
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for specific C-N cross-coupling reactions
                for name in C_N_COUPLING_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        has_aromatic_nucleophilic_substitution = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        # No break here, as multiple reactions might match and we want to record all

            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases for child (reaction)
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth, max_depth)

    # Start traversal. Assuming max_depth can be pre-calculated or is not needed.
    dfs_traverse(route, 1, 0)

    # Add structural constraint if the condition is met
    if has_aromatic_nucleophilic_substitution:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "C-N_cross-coupling",
                "operator": ">=",
                "value": 1
            }
        })

    return has_aromatic_nucleophilic_substitution, findings_json
