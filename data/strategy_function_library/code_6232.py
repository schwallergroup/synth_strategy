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


COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Buchwald-Hartwig",
    "Heck terminal vinyl",
    "Heck_terminal_vinyl",
    "Heck_non-terminal_vinyl",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Stille",
    "Negishi coupling",
    "Negishi",
    "Ullmann condensation",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "N-arylation_heterocycles",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final synthetic step is a cross-coupling reaction from a predefined list.
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

    has_late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_coupling, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                if len(reactants) >= 2:
                    # Check if this is a known coupling reaction
                    for rxn_type in COUPLING_REACTIONS:
                        if checker.check_reaction(rxn_type, rsmi):
                            has_late_stage_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            # Add structural constraint if found at the correct depth
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "cross-coupling reaction",
                                    "position": "last_or_penultimate_stage"
                                }
                            })
                            break
            except Exception:
                # Silently ignore nodes that don't conform to the expected reaction structure
                pass

        # Traverse children
        for child in node.get("children", []):
            # Optimization: stop traversing if we've already found the feature
            if has_late_stage_coupling:
                return
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_coupling, findings_json
