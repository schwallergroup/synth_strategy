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


HETEROCYCLES_OF_INTEREST = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "trioxane", "dioxepane", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "pyrrolidine", "piperidine", "piperazine", "morpholine",
    "thiomorpholine", "aziridine", "azetidine", "azepane", "diazepane",
    "indole", "quinoline", "isoquinoline", "purine", "carbazole", "acridine",
    "thiophene", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dithiane", "dithiolane", "benzothiophene", "oxathiolane", "dioxathiolane",
    "thiazolidine", "oxazolidine", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "benzimidazole", "benzoxazole", "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves early-stage formation (at a depth >= 2) of a specific heterocycle. The function checks for the appearance of a heterocycle from a predefined list (HETEROCYCLES_OF_INTEREST) in the product, ensuring it was not present in any of the reactants.
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

    early_heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal early_heterocycle_formation, findings_json

        if early_heterocycle_formation:
            return

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Only consider early stages (depth >= 2)
                if depth >= 2:
                    # Check for heterocycles in the product that are on our list
                    product_heterocycles = [
                        h for h in HETEROCYCLES_OF_INTEREST
                        if checker.check_ring(h, product_smiles)
                    ]

                    if product_heterocycles:
                        # For each new heterocycle, confirm it wasn't in the reactants
                        for heterocycle in product_heterocycles:
                            is_in_reactants = any(
                                checker.check_ring(heterocycle, r) for r in reactants_smiles
                            )
                            if not is_in_reactants:
                                early_heterocycle_formation = True
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                # Add the structural constraint if the condition is met
                                findings_json["structural_constraints"].append({
                                    "type": "positional",
                                    "details": {
                                        "target": "ring_formation",
                                        "position": "early_stage (depth >= 2)"
                                    }
                                })
                                return  # Strategy detected, stop this traversal path
            except (KeyError, IndexError):
                # Handle cases where rsmi might be missing or malformed
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return early_heterocycle_formation, findings_json
