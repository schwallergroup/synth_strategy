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
    Detects a strategy where a pyrimidine ring is formed while preserving a nitrile group. The strategy is flagged if the synthesis contains a reaction where a pyrimidine ring is newly formed on a reactant that already possesses a nitrile group, and this nitrile group is present in at least two distinct molecules across the entire synthetic route.
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

    # Track molecules with nitrile and pyrimidine
    molecules_with_nitrile = set()
    molecules_with_pyrimidine = set()

    # Track if we found a reaction that forms pyrimidine while preserving nitrile
    pyrimidine_formation_with_nitrile_preservation = False

    # Find the target molecule (final product)
    target_smiles = None
    if route["type"] == "mol":
        target_smiles = route["smiles"]

    def dfs_traverse(node, depth=0):
        nonlocal pyrimidine_formation_with_nitrile_preservation, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for nitrile and pyrimidine in this molecule
            has_nitrile = checker.check_fg("Nitrile", mol_smiles)
            has_pyrimidine = checker.check_ring("pyrimidine", mol_smiles)

            if has_nitrile:
                molecules_with_nitrile.add(mol_smiles)
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                print(f"Found molecule with nitrile at depth {depth}: {mol_smiles}")

            if has_pyrimidine:
                molecules_with_pyrimidine.add(mol_smiles)
                if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                print(f"Found molecule with pyrimidine at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has pyrimidine and nitrile
                product_has_pyrimidine = checker.check_ring("pyrimidine", product_smiles)
                product_has_nitrile = checker.check_fg("Nitrile", product_smiles)

                # Check reactants
                reactants_with_pyrimidine = []
                reactants_with_nitrile = []

                for r_smiles in reactants_smiles:
                    if checker.check_ring("pyrimidine", r_smiles):
                        reactants_with_pyrimidine.append(r_smiles)
                    if checker.check_fg("Nitrile", r_smiles):
                        reactants_with_nitrile.append(r_smiles)

                # In forward synthesis: reactants don't have pyrimidine, product does = pyrimidine formation
                if not any(reactants_with_pyrimidine) and product_has_pyrimidine:
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    print(f"Pyrimidine formation detected at depth {depth}")

                    # Check if nitrile is preserved through this reaction
                    if product_has_nitrile and len(reactants_with_nitrile) > 0:
                        print(f"Nitrile preserved during pyrimidine formation at depth {depth}")
                        pyrimidine_formation_with_nitrile_preservation = True
                        # Record the specific structural constraint for pyrimidine formation with nitrile preservation
                        if {"type": "co-occurrence", "details": {"description": "A pyrimidine ring formation reaction must occur where a Nitrile group is present in at least one reactant and is preserved in the product.", "targets": ["ring_formation", "pyrimidine", "Nitrile"]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "A pyrimidine ring formation reaction must occur where a Nitrile group is present in at least one reactant and is preserved in the product.", "targets": ["ring_formation", "pyrimidine", "Nitrile"]}})

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not a reaction (i.e., it's a chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if the target molecule has both nitrile and pyrimidine
    target_has_both = False
    target_has_pyrimidine = False
    target_has_nitrile = False

    if target_smiles:
        target_has_pyrimidine = checker.check_ring("pyrimidine", target_smiles)
        target_has_nitrile = checker.check_fg("Nitrile", target_smiles)
        target_has_both = target_has_pyrimidine and target_has_nitrile

        if target_has_both:
            print(f"Target molecule has both nitrile and pyrimidine: {target_smiles}")
            if {"type": "co-occurrence", "details": {"targets": ["pyrimidine", "Nitrile"], "position": "last_stage"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["pyrimidine", "Nitrile"], "position": "last_stage"}})
        elif target_has_pyrimidine:
            print(f"Target molecule has pyrimidine but no nitrile: {target_smiles}")
        elif target_has_nitrile:
            print(f"Target molecule has nitrile but no pyrimidine: {target_smiles}")

    # Strategy is present if:
    # 1. We found a reaction that forms pyrimidine while preserving nitrile
    # 2. Nitrile is present in multiple molecules throughout the synthesis
    # 3. Either the target has both features, or both features appear in the synthesis
    strategy_present = (
        pyrimidine_formation_with_nitrile_preservation
        and len(molecules_with_nitrile) >= 2
        and (
            target_has_both
            or (len(molecules_with_pyrimidine) > 0 and len(molecules_with_nitrile) > 0)
        )
    )

    if len(molecules_with_nitrile) >= 2:
        if {"type": "count", "details": {"target": "molecules_with_Nitrile", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "molecules_with_Nitrile", "operator": ">=", "value": 2}})

    if len(molecules_with_pyrimidine) > 0 and len(molecules_with_nitrile) > 0:
        if {"type": "co-occurrence", "details": {"description": "The route must contain at least one molecule with a pyrimidine ring and at least one molecule with a Nitrile group.", "targets": ["pyrimidine", "Nitrile"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "The route must contain at least one molecule with a pyrimidine ring and at least one molecule with a Nitrile group.", "targets": ["pyrimidine", "Nitrile"]}})

    if strategy_present:
        print("Detected nitrile-preserved pyrimidine formation strategy")
    else:
        print("Did not detect nitrile-preserved pyrimidine formation strategy")
        print(
            f"Pyrimidine formation with nitrile preservation: {pyrimidine_formation_with_nitrile_preservation}"
        )
        print(f"Number of molecules with nitrile: {len(molecules_with_nitrile)}")
        print(f"Target has both nitrile and pyrimidine: {target_has_both}")

    return strategy_present, findings_json
