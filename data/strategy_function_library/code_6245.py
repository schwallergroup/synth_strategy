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


COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "{Suzuki}",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Stille reaction_benzyl",
    "Stille reaction_allyl",
    "Stille reaction_aryl OTf",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Heck terminal vinyl",
    "Heck_terminal_vinyl",
    "Heck_non-terminal_vinyl",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
    "Ullmann condensation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves a late-stage coupling of two complex fragments,
    where each fragment contains multiple rings.
    """
    print("Analyzing route for late-stage fragment coupling")
    has_late_stage_complex_coupling = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_complex_coupling, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Late-stage includes final reaction
            print(f"Analyzing potential late-stage reaction at depth {depth}")
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Reaction SMILES: {rsmi}")

            # Check if this is a coupling reaction
            is_coupling = False
            detected_reaction = None

            # First try direct reaction checking
            for reaction_type in COUPLING_REACTIONS_OF_INTEREST:
                if checker.check_reaction(reaction_type, rsmi):
                    is_coupling = True
                    detected_reaction = reaction_type
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    print(f"Detected coupling reaction: {reaction_type}")
                    break

            if not is_coupling:
                print("Not a coupling reaction, skipping")
                return

            # If it's a coupling reaction and late-stage, add the constraint
            if depth <= 1:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "coupling_reaction",
                        "position": "late_stage",
                        "max_depth_from_product": 1
                    }
                })

            # Get reactants
            reactants = [r for r in rsmi.split(">")[0].split(".") if r and r.strip()]
            print(f"Found {len(reactants)} reactants")

            # Count complex fragments (fragments with 2+ rings)
            complex_fragments = 0
            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol:
                    # Count rings
                    ring_info = mol.GetRingInfo()
                    num_rings = len(ring_info.AtomRings())
                    print(f"Reactant {r} has {num_rings} rings")

                    if num_rings >= 2:  # Consider a fragment complex if it has 2+ rings
                        complex_fragments += 1
                        print(f"Found complex fragment with {num_rings} rings: {r}")

            print(f"Total complex fragments: {complex_fragments}")
            if complex_fragments >= 2:
                has_late_stage_complex_coupling = True
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "reactants_with_2_or_more_rings",
                        "operator": ">=",
                        "value": 2,
                        "scope": "per_reaction"
                    }
                })
                print(
                    f"Detected late-stage coupling of {complex_fragments} complex fragments at depth {depth}"
                )

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node['type'] != 'reaction': # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Result: {has_late_stage_complex_coupling}")
    return has_late_stage_complex_coupling, findings_json
