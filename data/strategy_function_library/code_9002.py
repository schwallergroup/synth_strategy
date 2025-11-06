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


import rdkit.Chem as Chem

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Heck",
    "Sonogashira",
    "Buchwald-Hartwig",
    "Ullmann-Goldberg",
    "N-arylation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage fragment coupling strategy. It first checks if the final step (depth=1) is a named reaction from a specific list: Suzuki, Negishi, Stille, Heck, Sonogashira, Buchwald-Hartwig, Ullmann-Goldberg, N-arylation, Hiyama-Denmark Coupling, Kumada cross-coupling, Aryllithium cross-coupling, decarboxylative_coupling.
    As a fallback, it identifies reactions where at least two structurally complex fragments (defined as having >12 atoms or >=2 rings) are joined.
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

    print("Starting late_stage_fragment_coupling_strategy analysis...")
    final_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal final_coupling_detected, findings_json

        if node["type"] == "reaction":
            # Extract reaction SMILES
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Check if this is the final synthetic step (depth=1)
                if depth == 1:
                    print(f"Analyzing late-stage reaction at depth {depth}: {rsmi}")

                    # Check for specific coupling reactions
                    for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Found late-stage coupling reaction: {rxn_type}")
                            final_coupling_detected = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "targets": [
                                        "Suzuki", "Negishi", "Stille", "Heck", "Sonogashira",
                                        "Buchwald-Hartwig", "Ullmann-Goldberg", "N-arylation",
                                        "Hiyama-Denmark Coupling", "Kumada cross-coupling",
                                        "Aryllithium cross-coupling", "decarboxylative_coupling"
                                    ],
                                    "position": "last_stage"
                                }
                            })
                            return

                    # If no specific coupling reaction found, check for complex fragments
                    if not final_coupling_detected and len(reactants) >= 2:
                        complex_reactants_count = 0
                        for reactant in reactants:
                            try:
                                mol = Chem.MolFromSmiles(reactant)
                                if mol:
                                    # Define complex fragment: >12 atoms or >=2 rings
                                    num_atoms = mol.GetNumAtoms()
                                    num_rings = mol.GetRingInfo().NumRings()

                                    if num_atoms > 12 or num_rings >= 2:
                                        complex_reactants_count += 1
                                        print(
                                            f"Found complex fragment: {reactant} (atoms: {num_atoms}, rings: {num_rings})"
                                        )
                            except Exception as e:
                                print(f"Error analyzing reactant {reactant}: {e}")

                        if complex_reactants_count >= 2:
                            print(
                                f"Found late-stage fragment coupling with {complex_reactants_count} complex fragments"
                            )
                            final_coupling_detected = True
                            findings_json["atomic_checks"]["named_reactions"].append("complex_fragment_coupling")
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "complex_fragment_coupling",
                                    "position": "last_stage"
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "complex_reactants",
                                    "operator": ">=",
                                    "value": 2
                                }
                            })
                            return

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage fragment coupling detected: {final_coupling_detected}")
    return final_coupling_detected, findings_json
