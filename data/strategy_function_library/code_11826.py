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


PYRAZOLE_FORMATION_REACTIONS = [
    "pyrazole",
    "{pyrazole}",
    "Pyrazole formation",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "[3+2]-cycloaddition of hydrazone and alkyne",
    "[3+2]-cycloaddition of hydrazone and alkene",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis strategy involves forming a pyrazole ring
    in the final step by checking against a specific list of formation reactions.
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

    final_product_has_pyrazole = False
    final_step_forms_pyrazole = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_pyrazole, final_step_forms_pyrazole, findings_json

        # Check if final product has pyrazole
        if node["type"] == "mol" and depth == 0:
            mol_smiles = node["smiles"]
            if checker.check_ring("pyrazole", mol_smiles):
                final_product_has_pyrazole = True
                findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "pyrazole",
                        "position": "final_product"
                    }
                })
                print(f"Final product contains pyrazole: {mol_smiles}")

        # Check if final reaction forms pyrazole
        elif node["type"] == "reaction" and depth <= 1:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product has pyrazole
            if checker.check_ring("pyrazole", product_smiles):
                print(f"Product has pyrazole: {product_smiles}")
                # This check is implicitly covered by the final_product_has_pyrazole check if it's the final product
                # but if it's an intermediate product of the last step, we still record the ring system.
                if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

                # Check if any reactant has pyrazole
                reactant_has_pyrazole = False
                for reactant_smiles in reactants_smiles:
                    if checker.check_ring("pyrazole", reactant_smiles):
                        reactant_has_pyrazole = True
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "pyrazole",
                                "scope": "reactants_of_last_stage"
                            }
                        })
                        print(f"Reactant already has pyrazole: {reactant_smiles}")
                        break

                # Check if this is a known pyrazole formation reaction
                is_pyrazole_reaction = False
                for name in PYRAZOLE_FORMATION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_pyrazole_reaction = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if not reactant_has_pyrazole and is_pyrazole_reaction:
                    final_step_forms_pyrazole = True
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "any_of_PYRAZOLE_FORMATION_REACTIONS",
                            "position": "last_stage"
                        }
                    })
                    print(f"Final step forms pyrazole through known reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            # If current node is a reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Final product has pyrazole: {final_product_has_pyrazole}")
    print(f"Final step forms pyrazole: {final_step_forms_pyrazole}")

    # Return True if pyrazole is formed in the final step
    result = final_product_has_pyrazole and final_step_forms_pyrazole
    return result, findings_json
