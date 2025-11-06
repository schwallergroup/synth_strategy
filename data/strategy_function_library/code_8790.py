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


ETHER_FORMATION_REACTIONS = [
    "Williamson Ether Synthesis",
    "{Williamson ether}",
    "Mitsunobu aryl ether",
    "Alcohol to ether",
    "Chan-Lam etherification",
    "Chan-Lam alcohol",
    "Ullmann-Goldberg Substitution aryl alcohol",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy combining pyrazole ring formation with a late-stage ether coupling, identified by specific named reactions.
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

    # Track if we found pyrazole formation and ether coupling
    pyrazole_formation_found = False
    late_stage_ether_coupling_found = False

    # Define a helper function to traverse the synthesis route
    def traverse_route(node, depth=0):
        nonlocal pyrazole_formation_found, late_stage_ether_coupling_found, findings_json

        # Check if this is a reaction node
        if node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]
            product_smiles = rxn_smiles.split(">")[-1]
            reactants_smiles = rxn_smiles.split(">")[0].split(".")

            # Check for pyrazole formation
            if not pyrazole_formation_found:
                # Check if product contains pyrazole ring and reactants do not
                if checker.check_ring("pyrazole", product_smiles):
                    reactants_have_pyrazole = any(
                        checker.check_ring("pyrazole", reactant) for reactant in reactants_smiles
                    )
                    if not reactants_have_pyrazole:
                        print(f"Found pyrazole formation reaction: {rxn_smiles}")
                        pyrazole_formation_found = True
                        if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            # Check for late-stage ether coupling
            # Late stage is defined as one of the final 3 synthetic steps.
            if depth <= 3:
                # Check for specific ether formation reactions by name
                for r in ETHER_FORMATION_REACTIONS:
                    if checker.check_reaction(r, rxn_smiles):
                        print(f"Found late-stage ether coupling reaction: {rxn_smiles}")
                        late_stage_ether_coupling_found = True
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break

        # Process molecule nodes
        elif node.get("type") == "mol":
            # Nothing specific to check in molecule nodes for this strategy
            pass

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical (mol), depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            traverse_route(child, new_depth)

    # Start traversal from the root
    traverse_route(route)

    combined_strategy = pyrazole_formation_found and late_stage_ether_coupling_found
    print(f"Combined pyrazole formation with ether coupling strategy: {combined_strategy}")
    print(f"Pyrazole formation found: {pyrazole_formation_found}")
    print(f"Late-stage ether coupling found: {late_stage_ether_coupling_found}")

    if pyrazole_formation_found and late_stage_ether_coupling_found:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "pyrazole_ring_formation",
                    "late_stage_ether_coupling"
                ]
            }
        })
    
    if late_stage_ether_coupling_found:
        # This constraint is met if any ether coupling is found within the specified depth.
        # The specific reactions are already added to atomic_checks if found.
        # We add the general positional constraint if the condition for late-stage is met.
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Williamson Ether Synthesis",
                    "{Williamson ether}",
                    "Mitsunobu aryl ether",
                    "Alcohol to ether",
                    "Chan-Lam etherification",
                    "Chan-Lam alcohol",
                    "Ullmann-Goldberg Substitution aryl alcohol",
                    "Ullmann condensation"
                ],
                "position": "within_final_4_stages"
            }
        })

    return combined_strategy, findings_json