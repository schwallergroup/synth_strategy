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


# Refactoring for Enumeration: List moved outside the function
HETEROCYCLE_TYPES_OF_INTEREST = [
    "pyrazole",
    "thiazole",
    "pyridine",
    "furan",
    "pyrrole",
    "imidazole",
    "oxazole",
    "isoxazole",
    "pyrimidine",
    "pyrazine",
    "triazole",
    "tetrazole",
    "indole",
    "benzimidazole",
    "benzothiazole",
    "benzoxazole",
    "quinoline",
    "isoquinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy characterized by the coupling of two or more fragments containing heterocycles from the HETEROCYCLE_TYPES_OF_INTEREST list. The overall route must also feature either the formation of a THP protecting group or a late-stage vinyl group introduction.
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

    # Initialize tracking variables
    has_convergent_step = False
    has_thp_protection = False
    has_late_vinyl_introduction = False
    has_heterocycle_coupling = False

    # List of heterocycles is now a module-level constant

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_step, has_thp_protection, has_late_vinyl_introduction, has_heterocycle_coupling, findings_json

        if node["type"] == "reaction":
            # Extract reactants and products
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for convergent synthesis (multiple complex reactants at depth > 2)
                if depth >= 2 and len(reactants_smiles) > 1:
                    # Check if at least two reactants are complex (contain heterocycles)
                    complex_reactants = 0
                    for reactant in reactants_smiles:
                        has_heterocycle = False
                        for heterocycle in HETEROCYCLE_TYPES_OF_INTEREST:
                            if checker.check_ring(heterocycle, reactant):
                                has_heterocycle = True
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                        if has_heterocycle:
                            complex_reactants += 1

                    if complex_reactants >= 2:
                        print(f"Found convergent step at depth {depth}")
                        has_convergent_step = True
                        if "convergent_heterocycle_step" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("convergent_heterocycle_step")

                # Check for THP protection
                if checker.check_ring("tetrahydropyran", product_smiles):
                    if "tetrahydropyran" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("tetrahydropyran")
                    # Check if this is a protection reaction (i.e., THP is formed)
                    if any(
                        not checker.check_ring("tetrahydropyran", reactant)
                        for reactant in reactants_smiles
                    ):
                        print(f"Found THP protection at depth {depth}")
                        has_thp_protection = True
                        if "thp_protection" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("thp_protection")

                # Check for late-stage vinyl introduction
                if depth <= 2:  # Late stage = low depth
                    product_has_vinyl = checker.check_fg("Vinyl", product_smiles)
                    if product_has_vinyl:
                        if "Vinyl" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Vinyl")

                        # Check if vinyl is being introduced (not present in all reactants)
                        reactants_without_vinyl = [
                            reactant
                            for reactant in reactants_smiles
                            if not checker.check_fg("Vinyl", reactant)
                        ]

                        if reactants_without_vinyl and len(reactants_without_vinyl) < len(
                            reactants_smiles
                        ):
                            print(f"Found vinyl group introduction at depth {depth}")
                            has_late_vinyl_introduction = True
                            if "vinyl_introduction" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("vinyl_introduction")

                # Check for heterocycle coupling at any depth
                product_heterocycle_count = 0
                for heterocycle in HETEROCYCLE_TYPES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_heterocycle_count += 1
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                if product_heterocycle_count >= 2:
                    # Check if this is a coupling reaction
                    reactant_heterocycle_counts = [
                        sum(
                            1
                            for heterocycle in HETEROCYCLE_TYPES_OF_INTEREST
                            if checker.check_ring(heterocycle, reactant)
                        )
                        for reactant in reactants_smiles
                    ]

                    # If product has more heterocycles than any single reactant, it's a coupling
                    if all(
                        count < product_heterocycle_count for count in reactant_heterocycle_counts
                    ):
                        print(f"Found heterocycle coupling at depth {depth}")
                        has_heterocycle_coupling = True
                        if "heterocycle_coupling" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("heterocycle_coupling")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if the strategy is detected
    result = (
        has_convergent_step
        and has_heterocycle_coupling
        and (has_late_vinyl_introduction or has_thp_protection)
    )

    if has_late_vinyl_introduction:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "vinyl_introduction",
                "position": "late_stage_depth_le_2"
            }
        })

    if has_convergent_step and has_heterocycle_coupling:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "convergent_heterocycle_step",
                    "heterocycle_coupling"
                ]
            }
        })

    if has_thp_protection or has_late_vinyl_introduction:
        targets = []
        if has_thp_protection:
            targets.append("thp_protection")
        if has_late_vinyl_introduction:
            targets.append("late_stage_vinyl_introduction")
        if targets:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": targets,
                    "operator": ">=",
                    "value": 1
                }
            })

    print(f"Final result: {result}")
    print(
        f"Convergent step: {has_convergent_step}, Heterocycle coupling: {has_heterocycle_coupling}"
    )
    print(
        f"Late vinyl introduction: {has_late_vinyl_introduction}, THP protection: {has_thp_protection}"
    )

    return result, findings_json
