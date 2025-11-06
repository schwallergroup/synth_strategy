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
    Detects a strategy involving two distinct nitrile formations:
    1. From halide substitution
    2. From amide conversion
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

    nitrile_from_halide = False
    nitrile_from_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_from_halide, nitrile_from_amide, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for nitrile formation from halide
                halide_types = ["Primary halide", "Secondary halide", "Tertiary halide", "Aromatic halide"]
                has_halide_reactant = False
                for h_type in halide_types:
                    if any(checker.check_fg(h_type, r) for r in reactants):
                        has_halide_reactant = True
                        if h_type not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(h_type)

                has_nitrile_product = checker.check_fg("Nitrile", product)
                if has_nitrile_product and "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                has_nitrile_reactant = any(checker.check_fg("Nitrile", r) for r in reactants)

                if has_halide_reactant and has_nitrile_product and not has_nitrile_reactant:
                    print(f"Found nitrile formation from halide at depth {depth}")
                    nitrile_from_halide = True
                    if "Nitrile formation from halide" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Nitrile formation from halide")

                # Check for nitrile formation from amide
                has_amide_reactant = any(
                    checker.check_fg("Primary amide", r)
                    for r in reactants
                )
                if has_amide_reactant and "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amide")

                if has_amide_reactant and has_nitrile_product and not has_nitrile_reactant:
                    print(f"Found nitrile formation from amide at depth {depth}")
                    nitrile_from_amide = True
                    if "Nitrile formation from amide" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Nitrile formation from amide")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Nitrile from halide: {nitrile_from_halide}, Nitrile from amide: {nitrile_from_amide}")

    result = nitrile_from_halide and nitrile_from_amide

    if result:
        # Add the structural constraint if both conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Nitrile formation from halide",
                    "Nitrile formation from amide"
                ]
            }
        })

    return result, findings_json