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
    "furan",
    "pyran",
    "dioxane",
    "tetrahydrofuran",
    "tetrahydropyran",
    "oxirane",
    "oxetane",
    "oxolane",
    "oxane",
    "dioxolane",
    "dioxolene",
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "pyrrolidine",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
    "indole",
    "quinoline",
    "isoquinoline",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "thiophene",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the late-stage (final or penultimate step) formation of a heterocycle from a predefined list (HETEROCYCLES_OF_INTEREST). The function confirms that a heterocycle from the list is present in the reaction product but absent from all reactants.
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

    final_steps = []
    heterocycle_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed

        if node["type"] == "reaction" and depth <= 1:
            final_steps.append((node, depth))

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    for step, depth in final_steps:
        if "rsmi" in step.get("metadata", {}):
            rsmi = step["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check which heterocycles are in the product
            product_heterocycles = []
            for heterocycle in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(heterocycle, product_smiles):
                    product_heterocycles.append(heterocycle)
                    # Record atomic check for ring system
                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

            if not product_heterocycles:
                continue

            # Check which heterocycles are in the reactants
            reactant_heterocycles = set()
            for reactant in reactants_smiles:
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, reactant):
                        reactant_heterocycles.add(heterocycle)

            # Check if any new heterocycles were formed
            new_heterocycles = [h for h in product_heterocycles if h not in reactant_heterocycles]

            if new_heterocycles:
                heterocycle_formed = True
                # Record structural constraint
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "ring_formation",
                        "position": "final or penultimate stage"
                    }
                })
                return True, findings_json

    return heterocycle_formed, findings_json