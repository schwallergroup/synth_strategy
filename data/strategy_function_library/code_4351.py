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


SUZUKI_REACTION_NAMES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects synthetic routes that feature a Suzuki coupling. It also identifies the presence of sulfonamide functional groups and tracks other reaction types. This function uses a predefined list of named Suzuki reaction variants for detection.
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

    # Track if we found Suzuki coupling
    found_suzuki = False
    suzuki_depth = None

    # Track if fluorine and sulfonamide are preserved
    fluorine_depths = set()
    sulfonamide_depths = set()

    # Track reaction types at each depth
    reaction_types = {}

    # Track biaryl formation
    biaryl_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki, suzuki_depth, biaryl_formed, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for sulfonamide
            if checker.check_fg("Sulfonamide", mol_smiles):
                sulfonamide_depths.add(depth)
                findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                print(f"Found sulfonamide at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for Suzuki coupling with expanded detection
                is_suzuki = False
                for name in SUZUKI_REACTION_NAMES:
                    if checker.check_reaction(name, rsmi):
                        is_suzuki = True
                        found_suzuki = True
                        suzuki_depth = depth
                        if depth not in reaction_types:
                            reaction_types[depth] = []
                        reaction_types[depth].append("suzuki")
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        print(f"Found Suzuki coupling at depth {depth}")
                        break

                if not is_suzuki:
                    # Check for amide formation
                    is_amide_formation = False
                    for name in AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(name, rsmi):
                            is_amide_formation = True
                            if depth not in reaction_types:
                                reaction_types[depth] = []
                            reaction_types[depth].append("amide_formation")
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                            print(f"Found amide formation at depth {depth}")
                            break

                    if not is_amide_formation:
                        # Check for esterification
                        if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                            if depth not in reaction_types:
                                reaction_types[depth] = []
                            reaction_types[depth].append("esterification")
                            findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")
                            print(f"Found esterification at depth {depth}")

                        # Check for any other reaction to track linear build-up
                        else:
                            if depth not in reaction_types:
                                reaction_types[depth] = []
                            reaction_types[depth].append("other")
                            print(f"Found other reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if this is a linear build-up strategy with Suzuki coupling
    # and preservation of functional groups

    # Linear build-up should have multiple reactions at different depths
    is_linear = len(reaction_types) >= 1
    if is_linear:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_reaction_step",
                "operator": ">=",
                "value": 1,
                "description": "The route must contain at least one reaction step."
            }
        })

    # Functional groups should be preserved throughout synthesis
    # Note: The original code's fluorine_preserved logic is based on `fluorine_depths` which is never populated.
    # This part of the structural constraint will only be added if the condition is met, which it currently won't be.
    fluorine_preserved = len(fluorine_depths) >= 2
    if fluorine_preserved:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "stages_with_fluorine",
                "operator": ">=",
                "value": 2,
                "description": "The route logic requires fluorine to be present in at least two different stages, although the implementation for detecting fluorine is missing in the source code."
            }
        })

    # Either sulfonamide is preserved or it's not part of the strategy
    sulfonamide_present = len(sulfonamide_depths) > 0
    sulfonamide_preserved = len(sulfonamide_depths) >= 1 if sulfonamide_present else True

    # For a true linear build-up with Suzuki strategy, we need:
    # 1. Suzuki coupling present
    # 2. Linear build-up (multiple reactions)
    # 3. Fluorine preservation
    # 4. Sulfonamide preservation (if present)
    # 5. Ideally, biaryl formation via Suzuki

    strategy_present = found_suzuki and is_linear and fluorine_preserved and sulfonamide_preserved

    if found_suzuki:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki"
                ],
                "min_count": 1,
                "description": "At least one of the specified Suzuki reactions must be present in the route."
            }
        })

    print(f"Linear build-up with Suzuki strategy detected: {strategy_present}")
    print(f"Reaction types by depth: {reaction_types}")
    print(f"Fluorine depths: {fluorine_depths}")
    print(f"Sulfonamide depths: {sulfonamide_depths}")
    print(f"Found Suzuki: {found_suzuki}")
    print(f"Is linear: {is_linear}")
    print(f"Fluorine preserved: {fluorine_preserved}")
    print(f"Sulfonamide preserved: {sulfonamide_preserved}")

    return strategy_present, findings_json
