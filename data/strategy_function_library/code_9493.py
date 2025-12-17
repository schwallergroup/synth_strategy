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


# Refactored lists at module level
BICYCLIC_HETEROCYCLES = [
    "indole",
    "quinoline",
    "isoquinoline",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "purine",
    "carbazole",
    "acridine",
    "benzothiophene",
    "pteridin",
    "indazole",
]

HETEROCYCLE_FORMING_REACTIONS = [
    "Pictet-Spengler",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "Niementowski_quinazoline",
    "tetrazole_terminal",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "3-nitrile-pyridine",
    "pyrazole",
    "Fischer indole",
    "Friedlaender chinoline",
    "benzofuran",
    "benzothiophene",
    "indole",
    "oxadiazole",
    "imidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the early-stage formation of specific bicyclic heterocyclic cores. This is determined by checking for the presence of a core from the `BICYCLIC_HETEROCYCLES` list in the product, or by identifying a reaction from the `HETEROCYCLE_FORMING_REACTIONS` list, and ensuring the core was not present in the reactants. The check is restricted to reactions occurring at a depth of 2 or greater.
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

    core_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal core_formation_detected, findings_json

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an early-stage reaction (depth >= 2)
            # In retrosynthetic analysis, higher depth means earlier stage
            if depth >= 2:
                print(f"Checking early-stage reaction at depth {depth}: {rsmi}")
                # Add positional constraint if depth >= 2
                if {"type": "positional", "details": {"target": "bicyclic_heterocycle_formation", "position": "early_stage (depth >= 2)"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "bicyclic_heterocycle_formation", "position": "early_stage (depth >= 2)"}})

                # Check for bicyclic heterocyclic system in product
                product_has_core = False
                detected_cores = []

                # Check for bicyclic heterocycles
                for ring in BICYCLIC_HETEROCYCLES:
                    if checker.check_ring(ring, product):
                        product_has_core = True
                        detected_cores.append(ring)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        print(f"Found heterocyclic core: {ring} in product")

                # Check for heterocycle-forming reactions
                reaction_forms_heterocycle = False
                for rxn_type in HETEROCYCLE_FORMING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected heterocycle-forming reaction: {rxn_type}")
                        reaction_forms_heterocycle = True
                        product_has_core = True # A reaction forming a heterocycle implies product has core
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if product_has_core or reaction_forms_heterocycle:
                    # Check if reactants don't have the same heterocyclic system
                    reactants_have_core = False
                    for reactant in reactants:
                        # Check for bicyclic cores in reactants
                        for ring in BICYCLIC_HETEROCYCLES:
                            if checker.check_ring(ring, reactant):
                                reactants_have_core = True
                                print(f"Reactant also contains the core: {ring} in {reactant}")
                                break

                        if reactants_have_core:
                            break

                    if not reactants_have_core:
                        print(f"Early heterocyclic core formation detected at depth {depth}")
                        core_formation_detected = True
                        # Add negation constraint if reactants don't have core
                        if {"type": "negation", "details": {"target": "bicyclic_heterocycle_in_reactants", "scope": "reaction_step_with_potential_formation"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "bicyclic_heterocycle_in_reactants", "scope": "reaction_step_with_potential_formation"}})
                else:
                    print(f"No heterocyclic core found in product: {product}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return core_formation_detected, findings_json
