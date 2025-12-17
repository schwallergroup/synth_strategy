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


LATE_STAGE_ETHER_COUPLINGS = [
    "Williamson Ether Synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Mitsunobu aryl ether",
    "Mitsunobu aryl ether (intramolecular)",
    "Chan-Lam etherification",
    "Ullmann-Goldberg Substitution aryl alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage ether bond formation between fragments, identified by a specific list of named reactions.
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

    late_ether_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_ether_coupling, findings_json
        if late_ether_coupling:  # Optimization: stop traversing if already found
            return

        # Consider reactions at depth 0 or 1 as late-stage
        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check if we have multiple reactants (fragment coupling)
                if len(reactants_smiles) >= 2:
                    # Check for specific ether formation reactions
                    for reaction_name in LATE_STAGE_ETHER_COUPLINGS:
                        if checker.check_reaction(reaction_name, rsmi):
                            late_ether_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            # Add structural constraints when conditions are met
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "ether_coupling_reaction",
                                    "position": "late_stage (depth <= 1)"
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "reactants_in_ether_coupling",
                                    "operator": ">=",
                                    "value": 2
                                }
                            })
                            return
            except Exception:
                # Silently ignore errors in reaction analysis
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_ether_coupling, findings_json
