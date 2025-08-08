#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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


def main(route):
    """
    This function detects a synthetic strategy involving both cross-coupling (C-C bond formation)
    and SNAr reactions (C-N bond formation) in the same route.
    """
    has_cross_coupling = False
    has_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal has_cross_coupling, has_snar

        indent = "  " * depth
        print(f"{indent}Examining node type: {node['type']}")

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"{indent}Checking reaction SMILES: {rsmi}")

                # Check for cross-coupling reactions using checker functions
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Stille reaction_vinyl", rsmi)
                    or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                    or checker.check_reaction("Aryllithium cross-coupling", rsmi)
                    or checker.check_reaction("decarboxylative_coupling", rsmi)
                    or checker.check_reaction("{Suzuki}", rsmi)
                ):
                    has_cross_coupling = True
                    print(f"{indent}Cross-coupling reaction detected: {rsmi}")

                # If no direct cross-coupling reaction detected, check for characteristic patterns
                if not has_cross_coupling:
                    # Extract reactants and product
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant has boronic acid/ester and another has halide
                        has_boronic = any(
                            checker.check_fg("Boronic acid", r)
                            or checker.check_fg("Boronic ester", r)
                            for r in reactants
                        )
                        has_halide = any(
                            checker.check_fg("Aromatic halide", r)
                            or checker.check_fg("Primary halide", r)
                            or checker.check_fg("Secondary halide", r)
                            or checker.check_fg("Tertiary halide", r)
                            for r in reactants
                        )

                        if has_boronic and has_halide:
                            has_cross_coupling = True
                            print(
                                f"{indent}Cross-coupling reaction inferred from functional groups: {rsmi}"
                            )
                    except Exception as e:
                        print(f"{indent}Error analyzing reactants: {e}")

                # Check for SNAr reactions using checker functions
                if (
                    checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                    or checker.check_reaction("Goldberg coupling", rsmi)
                    or checker.check_reaction("Goldberg coupling aryl amine-aryl chloride", rsmi)
                    or checker.check_reaction("Goldberg coupling aryl amide-aryl chloride", rsmi)
                    or checker.check_reaction("Ullmann condensation", rsmi)
                    or checker.check_reaction("{Buchwald-Hartwig}", rsmi)
                    or checker.check_reaction("{N-arylation_heterocycles}", rsmi)
                ):
                    has_snar = True
                    print(f"{indent}SNAr reaction detected: {rsmi}")

                # If no direct SNAr reaction detected, check for characteristic patterns
                if not has_snar:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant has amine and another has aromatic halide
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Aniline", r)
                            for r in reactants
                        )
                        has_aromatic_halide = any(
                            checker.check_fg("Aromatic halide", r) for r in reactants
                        )

                        # Check if product has C-N bond in aromatic system
                        if has_amine and has_aromatic_halide:
                            has_snar = True
                            print(f"{indent}SNAr reaction inferred from functional groups: {rsmi}")
                    except Exception as e:
                        print(f"{indent}Error analyzing reactants: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting route traversal...")
    dfs_traverse(route)

    # Debug output
    print(f"Cross-coupling detected: {has_cross_coupling}")
    print(f"SNAr detected: {has_snar}")

    # Return True if both cross-coupling and SNAr are detected
    return has_cross_coupling and has_snar
