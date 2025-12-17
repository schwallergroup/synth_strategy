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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects the formation of a cyclopropane ring from an alkene precursor during the middle stages of a synthesis. The 'middle stage' is defined as any step that is not one of the first two steps of the overall synthesis.
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

    found_cyclopropanation = False

    # Calculate maximum depth of the route
    def calculate_max_depth(node, current_depth=0):
        if not node.get("children", []):
            return current_depth

        max_child_depth = 0
        for child in node.get("children", []):
            child_depth = calculate_max_depth(child, current_depth + 1)
            max_child_depth = max(max_child_depth, child_depth)

        return max_child_depth

    max_depth = calculate_max_depth(route)

    # Define what constitutes "middle" of synthesis
    mid_start = 1  # Excludes step 0 (the root molecule)
    mid_end = max(3, max_depth - 2)  # Excludes the first two synthesis steps in longer routes

    def dfs_traverse(node, depth=0):
        nonlocal found_cyclopropanation, findings_json

        if node["type"] == "reaction" and mid_start <= depth <= mid_end:  # Mid-synthesis
            # Record positional constraint if met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "ring_formation",
                    "position": "middle_stages"
                }
            })

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains cyclopropane
                cyclopropane_in_product = checker.check_ring("cyclopropane", product)
                if cyclopropane_in_product:
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")

                # Check if cyclopropane is newly formed (not in reactants)
                cyclopropane_in_reactants = False
                for reactant in reactants:
                    if checker.check_ring("cyclopropane", reactant):
                        cyclopropane_in_reactants = True
                        break

                # Check for alkene in reactants (typical substrate for cyclopropanation)
                alkene_in_reactants = False
                for reactant in reactants:
                    if checker.check_fg("Alkene", reactant):
                        alkene_in_reactants = True
                        findings_json["atomic_checks"]["functional_groups"].append("Alkene")
                        break

                # Verify this is a cyclopropanation reaction
                if not cyclopropane_in_reactants: # Negation constraint met
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "cyclopropane",
                            "scope": "reactants_of_formation_step"
                        }
                    })
                    # The most robust check is the formation of a cyclopropane from an alkene.
                    if alkene_in_reactants and cyclopropane_in_product: # Co-occurrence constraint met
                        found_cyclopropanation = True
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "reaction_step",
                                "targets": [
                                    "ring_formation",
                                    "Alkene"
                                ]
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_cyclopropanation, findings_json
