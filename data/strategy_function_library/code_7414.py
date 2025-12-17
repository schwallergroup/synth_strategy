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


HYDRAZINE_DERIVATIVES = [
    "Hydrazine",
    "Hydrazone",
    "Acylhydrazine",
    "Hydrazone amide",
]
NN_HETEROCYCLES = [
    "pyrazole",
    "indazole",
    "triazole",
    "tetrazole",
    "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """ 
    This function is a placeholder for the updated description.
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

    found_hydrazine_cyclization = False

    def dfs_traverse(node, depth, max_depth):
        nonlocal found_hydrazine_cyclization, findings_json

        if found_hydrazine_cyclization:
            return  # Early return if already found

        if node["type"] == "reaction":
            reaction = node
            rsmi = reaction.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for specific named reactions known to be hydrazine-based heterocyclizations
            specific_reactions_to_check = [
                "pyrazole",
                "[3+2]-cycloaddition of hydrazone and alkyne",
                "[3+2]-cycloaddition of hydrazone and alkene",
                "Michael-induced ring closure from hydrazone",
                "[3+2]-cycloaddition of diazoalkane and alkyne",
                "[3+2]-cycloaddition of diazoalkane and alkene",
                "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
                "[3+2]-cycloaddition of diazoalkane and alpha-alkene"
            ]

            for r_name in specific_reactions_to_check:
                if checker.check_reaction(r_name, rsmi):
                    found_hydrazine_cyclization = True
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    return

            # Check for a hydrazine derivative in reactants leading to an N-N heterocycle in the product
            has_hydrazine_derivative = False
            found_fgs = []
            for fg in HYDRAZINE_DERIVATIVES:
                for r in reactants_smiles:
                    if checker.check_fg(fg, r):
                        has_hydrazine_derivative = True
                        if fg not in found_fgs:
                            found_fgs.append(fg)
                        break
            
            if has_hydrazine_derivative:
                for fg_name in found_fgs:
                    if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg_name)

                product_forms_nn_heterocycle = False
                found_rings = []
                for ring in NN_HETEROCYCLES:
                    if checker.check_ring(ring, product_smiles):
                        product_forms_nn_heterocycle = True
                        if ring not in found_rings:
                            found_rings.append(ring)
                        break

                if product_forms_nn_heterocycle:
                    found_hydrazine_cyclization = True
                    for ring_name in found_rings:
                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                    
                    # Add structural constraint if both conditions are met
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "hydrazine_derivative_reactant",
                                "nn_heterocycle_product"
                            ],
                            "scope": "single_reaction"
                        }
                    })
                    return

        # Traverse children, propagating context
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth, max_depth) # depth propagation is illustrative

    # Start traversal
    dfs_traverse(route, 1, 0) # max_depth calculation is assumed to be handled by a wrapper

    return found_hydrazine_cyclization, findings_json
