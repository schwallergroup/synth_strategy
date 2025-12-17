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


HETEROCYCLES_FROM_NITRILES = [
    "pyrazole",
    "tetrazole",
    "triazole",
    "oxadiazole",
    "thiadiazole",
    "isoxazole",
    "isothiazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyridazine",
    "pyrazine",
    "pyridine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a two-step pattern: 1) the presence of a nitrile-containing intermediate, and 2) its subsequent conversion into a specific heterocycle in a late-stage step (depth <= 3). The list of target heterocycles is defined in HETEROCYCLES_FROM_NITRILES.
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

    has_nitrile_intermediate = False
    nitrile_to_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitrile_intermediate, nitrile_to_heterocycle, findings_json

        if node["type"] == "mol":
            # Check if molecule is a nitrile intermediate
            if checker.check_fg("Nitrile", node["smiles"]) and depth > 0:
                has_nitrile_intermediate = True
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                print(f"Nitrile intermediate found at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check if reaction converts nitrile to heterocycle
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has nitrile
            nitrile_reactants = [r for r in reactants if r and checker.check_fg("Nitrile", r)]
            reactant_has_nitrile = len(nitrile_reactants) > 0
            if reactant_has_nitrile:
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

            # Check if product has any heterocyclic ring
            product_heterocycles_found = []
            for ring in HETEROCYCLES_FROM_NITRILES:
                if checker.check_ring(ring, product):
                    product_heterocycles_found.append(ring)
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
            product_has_heterocycle = len(product_heterocycles_found) > 0

            # Check if reactants already have the heterocycle
            reactants_combined = ".".join(reactants)
            reactants_have_heterocycle = any(
                checker.check_ring(ring, reactants_combined) for ring in HETEROCYCLES_FROM_NITRILES
            )

            # Check if this is a late-stage reaction (depth <= 3) that forms a heterocycle from a nitrile
            if (
                reactant_has_nitrile
                and product_has_heterocycle
                and not reactants_have_heterocycle
                and depth <= 3
            ):
                nitrile_to_heterocycle = True
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                
                # Add structural constraints if this condition is met
                if {"type": "co-occurrence", "details": {"targets": ["Nitrile", "ring_formation"]}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Nitrile", "ring_formation"]}})
                if {"type": "positional", "details": {"target": "Nitrile", "position": "not_last_stage"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitrile", "position": "not_last_stage"}})
                if {"type": "count", "details": {"target": "depth_of_ring_formation", "operator": "<=", "value": 3}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "depth_of_ring_formation", "operator": "<=", "value": 3}})

                print(f"Nitrile to heterocycle transformation found at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                print(f"Nitrile reactant: {nitrile_reactants[0] if nitrile_reactants else 'None'}")
                print(f"Heterocycle product: {product}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = has_nitrile_intermediate and nitrile_to_heterocycle
    print(f"Has nitrile intermediate: {has_nitrile_intermediate}")
    print(f"Has nitrile to heterocycle transformation: {nitrile_to_heterocycle}")
    print(f"Final result: {result}")

    return result, findings_json
