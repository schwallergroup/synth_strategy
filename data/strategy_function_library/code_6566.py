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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects convergent synthesis involving coupling of a heterocyclic core
    (specifically purine) with another fragment in a late-stage reaction.
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

    late_stage_coupling = False
    purine_involved = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling, purine_involved, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Late-stage reaction (depth 0, 1, or 2)
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have multiple fragments being combined
                if len(reactants) >= 2:
                    # Check for purine scaffold in reactants
                    purine_in_reactants = False
                    for reactant in reactants:
                        if checker.check_ring("purine", reactant):
                            purine_in_reactants = True
                            if "purine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("purine")
                            break

                    # Check if purine is preserved in product
                    purine_in_product = checker.check_ring("purine", product)
                    if purine_in_product and "purine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("purine")

                    # Only set purine_involved if purine is in both reactants and product
                    if purine_in_reactants and purine_in_product:
                        purine_involved = True

                        # Check if this is a coupling reaction
                        is_coupling = False

                        # Check for common coupling reaction types
                        if checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        ):
                            is_coupling = True
                            if "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)")

                        if is_coupling:
                            late_stage_coupling = True
                            # Add positional constraint if late_stage_coupling is true
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                                    "position": "late_stage"
                                }
                            })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = late_stage_coupling and purine_involved

    if result:
        # Add co-occurrence constraint if both conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "purine"
                ]
            }
        })

    return result, findings_json
