from typing import Tuple, Dict, List
import copy
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
    This function detects a synthetic strategy where an ether is formed via
    mesylate activation of an alcohol followed by nucleophilic substitution.
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

    # Track mesylate intermediates and ether formations
    mesylate_intermediates = set()
    ether_formations = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json, mesylate_intermediates, ether_formations
        # For reaction nodes, check for mesylate formation and ether formation
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for mesylate formation
            if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                print(f"Detected mesylate formation reaction at depth {depth}: {rsmi}")
                findings_json["atomic_checks"]["named_reactions"].append("Formation of Sulfonic Esters")
                # The product should contain a mesylate group
                if checker.check_fg("Mesylate", product):
                    mesylate_intermediates.add(product)
                    findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
            # Else, check if this is the ether formation step from a mesylate
            else:
                for reactant in reactants:
                    # If a reactant has a mesylate group
                    if checker.check_fg("Mesylate", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
                        # And the product has an ether group
                        if checker.check_fg("Ether", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")
                            # And the mesylate is not in the product (it was displaced)
                            if not checker.check_fg("Mesylate", product):
                                print(
                                    f"Detected unlabeled ether formation from mesylate at depth {depth}"
                                )
                                print(f"Mesylate reactant: {reactant}")
                                print(f"Ether product: {product}")
                                ether_formations.append((reactant, product))
                                # Found the reaction, no need to check other reactants
                                break

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both mesylate formation and ether formation using mesylate
    has_mesylate_formation = len(mesylate_intermediates) > 0
    has_ether_formation_via_mesylate = len(ether_formations) > 0

    print(f"Has mesylate formation: {has_mesylate_formation}")
    print(f"Has ether formation via mesylate: {has_ether_formation_via_mesylate}")
    print(f"Number of mesylate intermediates: {len(mesylate_intermediates)}")
    print(f"Number of ether formations via mesylate: {len(ether_formations)}")

    result = has_mesylate_formation and has_ether_formation_via_mesylate

    if result:
        # Add the structural constraint if both conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Formation of Sulfonic Esters",
                    "Mesylate",
                    "Ether"
                ],
                "description": "The route must contain both a 'Formation of Sulfonic Esters' reaction and a subsequent reaction where a 'Mesylate' is consumed to form an 'Ether'."
            }
        })

    return result, findings_json
