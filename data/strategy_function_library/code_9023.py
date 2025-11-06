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
    "Michael-induced ring closure from hydrazone",
    "[3+2]-cycloaddition of hydrazone and alkyne",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves late-stage formation of a pyrazole ring. This is confirmed if the final product contains a pyrazole ring that is formed in the final synthetic step. The formation reaction is identified from a predefined list of named reactions: 'pyrazole', 'Michael-induced ring closure from hydrazone', and '[3+2]-cycloaddition of hydrazone and alkyne'.
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
    first_reaction_forms_pyrazole = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_pyrazole, first_reaction_forms_pyrazole, findings_json

        if node["type"] == "mol" and depth == 0:
            # Check if final product has pyrazole
            if checker.check_ring("pyrazole", node["smiles"]):
                final_product_has_pyrazole = True
                findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "pyrazole",
                        "position": "final_product"
                    }
                })
                print(f"Final product contains pyrazole ring: {node['smiles']}")

        elif node["type"] == "reaction" and depth == 1:  # First reaction in retrosynthetic analysis
            # Check if this reaction forms pyrazole
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if the reaction is a pyrazole formation reaction
                is_pyrazole_reaction = False
                for r in PYRAZOLE_FORMATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        is_pyrazole_reaction = True
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break

                if is_pyrazole_reaction:
                    print("Reaction identified as potential pyrazole formation")

                    # Verify product has pyrazole but reactants don't
                    product_has_pyrazole = checker.check_ring("pyrazole", product)
                    reactants_have_pyrazole = any(
                        checker.check_ring("pyrazole", r) for r in reactants if r
                    )

                    if product_has_pyrazole and not reactants_have_pyrazole:
                        print("Confirmed: Product has pyrazole but reactants don't")
                        first_reaction_forms_pyrazole = True
                        if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                        
                        # Add the specific positional constraint for the reaction
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": [
                                    "pyrazole",
                                    "Michael-induced ring closure from hydrazone",
                                    "[3+2]-cycloaddition of hydrazone and alkyne"
                                ],
                                "position": "last_stage"
                            }
                        })
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "ring_formation",
                                "event_details": {
                                    "ring_name": "pyrazole"
                                },
                                "position": "last_stage"
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node['type'] == 'mol'
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    result = final_product_has_pyrazole and first_reaction_forms_pyrazole
    print(
        f"Final result: {result} (final_product_has_pyrazole={final_product_has_pyrazole}, first_reaction_forms_pyrazole={first_reaction_forms_pyrazole})"
    )
    return result, findings_json
