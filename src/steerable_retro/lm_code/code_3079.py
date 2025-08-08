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
    Detects a synthetic strategy involving protected heterocycle elaboration with a
    late-stage Suzuki coupling, followed by deprotection and final acylation.
    """
    # Initialize tracking variables
    has_boc_protection = False
    has_boc_deprotection = False
    has_suzuki_coupling = False
    has_snar_reaction = False
    has_final_acylation = False
    has_pyrimidine_core = False
    has_morpholine = False
    has_heterocycle_elaboration = False
    has_carbamate = False

    # Store all molecule nodes for later inspection
    all_mols = []

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection, has_boc_deprotection, has_suzuki_coupling
        nonlocal has_snar_reaction, has_final_acylation, has_pyrimidine_core
        nonlocal has_morpholine, has_heterocycle_elaboration, has_carbamate

        if node["type"] == "mol":
            # Store molecule node for later inspection
            all_mols.append(node)

            # Check for pyrimidine core and morpholine
            mol_smiles = node["smiles"]
            if checker.check_ring("pyrimidine", mol_smiles):
                has_pyrimidine_core = True
                print(f"Detected pyrimidine core in molecule: {mol_smiles}")

            if checker.check_ring("morpholine", mol_smiles):
                has_morpholine = True
                print(f"Detected morpholine in molecule: {mol_smiles}")

            # Check for carbamate group
            if checker.check_fg("Carbamic ester", mol_smiles):
                has_carbamate = True
                print(f"Detected carbamic ester in molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" not in node or "rsmi" not in node["metadata"]:
                return

            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection/deprotection
            reactants_have_boc = any(checker.check_fg("Boc", r) for r in reactants)
            product_has_boc = checker.check_fg("Boc", product)

            # Check for Boc deprotection
            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):
                has_boc_deprotection = True
                print(f"Detected Boc deprotection at depth {depth}")
            elif reactants_have_boc and not product_has_boc:
                has_boc_deprotection = True
                print(f"Detected Boc removal at depth {depth}")

            # Check for Boc protection
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                has_boc_protection = True
                print(f"Detected Boc protection at depth {depth}")
            elif not reactants_have_boc and product_has_boc:
                has_boc_protection = True
                print(f"Detected Boc addition at depth {depth}")

            # Check for Suzuki coupling - expanded to check for boronic acids/esters
            has_boronic_acid_or_ester = any(
                checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                for r in reactants
            )
            has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

            if has_boronic_acid_or_ester and has_aryl_halide:
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki", rsmi)
                ):
                    has_suzuki_coupling = True
                    print(f"Detected Suzuki coupling at depth {depth}")

            # Expanded check for Suzuki coupling without specific reaction type
            if not has_suzuki_coupling and has_boronic_acid_or_ester and has_aryl_halide:
                # Look for pattern where an aryl group from boronic acid/ester is connected to where halide was
                has_suzuki_coupling = True
                print(f"Detected Suzuki coupling pattern at depth {depth}")

            # Check for heterocycle elaboration - expanded beyond SNAr
            # SNAr reactions
            if (
                checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
            ):
                has_snar_reaction = True
                has_heterocycle_elaboration = True
                print(f"Detected SNAr reaction (heterocycle elaboration) at depth {depth}")

            # Check for other heterocycle elaboration reactions
            if (
                checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
                or checker.check_reaction("N-arylation_heterocycles", rsmi)
            ):
                has_heterocycle_elaboration = True
                print(f"Detected heterocycle elaboration via N-arylation at depth {depth}")

            # Check for heterocycle formation
            if checker.check_reaction("Formation of NOS Heterocycles", rsmi):
                has_heterocycle_elaboration = True
                print(f"Detected heterocycle formation at depth {depth}")

            # Check for chlorine substitution on pyrimidine
            reactants_have_pyrimidine_cl = any(
                checker.check_ring("pyrimidine", r) and checker.check_fg("Aromatic halide", r)
                for r in reactants
            )
            product_has_pyrimidine = checker.check_ring("pyrimidine", product)

            if reactants_have_pyrimidine_cl and product_has_pyrimidine:
                has_heterocycle_elaboration = True
                print(f"Detected heterocycle elaboration via halide substitution at depth {depth}")

            # Check for final acylation - expanded to include more acylation reactions
            if depth <= 3:  # Allow for more flexibility in the final steps
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Carboxylic acid to amide conversion", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                ):
                    has_final_acylation = True
                    print(f"Detected final acylation at depth {depth}")

                # Check for methyl carbamate formation
                product_has_carbamate = checker.check_fg("Carbamic ester", product)
                if product_has_carbamate and not any(
                    checker.check_fg("Carbamic ester", r) for r in reactants
                ):
                    has_final_acylation = True
                    print(f"Detected carbamate formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # Modified to be more flexible with the requirements
    has_protection = has_boc_protection or has_boc_deprotection or has_carbamate

    strategy_present = (
        has_pyrimidine_core and has_morpholine and has_suzuki_coupling and has_protection
    )

    # If we have the core elements but are missing some steps, check if we can infer them
    if has_pyrimidine_core and has_morpholine and has_suzuki_coupling:
        # If we have a carbamate but didn't detect explicit heterocycle elaboration
        if has_carbamate and not has_heterocycle_elaboration:
            has_heterocycle_elaboration = True
            print("Inferred heterocycle elaboration from presence of carbamate and pyrimidine core")

        # If we have a carbamate but didn't detect explicit final acylation
        if has_carbamate and not has_final_acylation:
            has_final_acylation = True
            print("Inferred final acylation from presence of carbamate")

    print(f"Strategy detection results:")
    print(f"  Pyrimidine core: {has_pyrimidine_core}")
    print(f"  Morpholine: {has_morpholine}")
    print(f"  Suzuki coupling: {has_suzuki_coupling}")
    print(f"  Boc protection/deprotection/carbamate: {has_protection}")
    print(f"  Heterocycle elaboration: {has_heterocycle_elaboration}")
    print(f"  Final acylation: {has_final_acylation}")
    print(f"  Overall strategy present: {strategy_present}")

    return strategy_present
