#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a synthetic strategy involving benzimidazole formation
    from a diamine and aldehyde, followed by late-stage cross-coupling to introduce
    structural diversity.
    """
    # Initialize flags for key transformations
    benzimidazole_formation = False
    late_stage_coupling = False
    nitro_reduction = False
    snar_reaction = False

    # Track molecules containing benzimidazole
    benzimidazole_molecules = set()

    def dfs_traverse(node, current_depth=0):
        nonlocal benzimidazole_formation, late_stage_coupling, nitro_reduction, snar_reaction

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[-1]

            # Get depth information
            depth_match = re.search(
                r"Depth: (\d+)",
                node.get("metadata", {}).get("ID", f"Depth: {current_depth}"),
            )
            depth = int(depth_match.group(1)) if depth_match else current_depth

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for benzimidazole formation
            if checker.check_ring("benzimidazole", product):
                print(f"Product contains benzimidazole at depth {depth}")
                benzimidazole_molecules.add(product)

                # Check if this is a benzimidazole formation reaction
                has_diamine = False
                has_aldehyde = False

                for reactant in reactants:
                    if checker.check_fg("Aldehyde", reactant):
                        has_aldehyde = True
                        print(f"Found aldehyde reactant: {reactant}")

                    # Check for diamine (o-phenylenediamine)
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        # Count primary amine groups
                        amine_count = 0
                        for atom in r_mol.GetAtoms():
                            if atom.GetSymbol() == "N" and atom.GetDegree() == 1:
                                amine_count += 1

                        if amine_count >= 2 and checker.check_fg(
                            "Primary amine", reactant
                        ):
                            has_diamine = True
                            print(f"Found diamine reactant: {reactant}")

                # Check if this is a benzimidazole formation reaction
                if has_aldehyde and has_diamine:
                    benzimidazole_formation = True
                    print(f"Detected benzimidazole formation at depth {depth}")

                # Alternative check using reaction type
                if checker.check_reaction("benzimidazole_derivatives_aldehyde", rsmi):
                    benzimidazole_formation = True
                    print(f"Detected benzimidazole formation reaction at depth {depth}")

            # Check for late-stage cross-coupling
            if depth <= 1:  # Late-stage coupling happens at low depths
                # Check if product contains benzimidazole
                if checker.check_ring("benzimidazole", product):
                    # Check for Suzuki coupling reaction
                    if checker.check_reaction("Suzuki", rsmi):
                        late_stage_coupling = True
                        print(f"Detected late-stage Suzuki coupling at depth {depth}")
                    elif checker.check_reaction(
                        "Suzuki coupling with boronic acids", rsmi
                    ):
                        late_stage_coupling = True
                        print(
                            f"Detected late-stage Suzuki coupling with boronic acids at depth {depth}"
                        )
                    else:
                        # Fallback to checking reactants for boronic acid and aromatic halide
                        has_boronic_acid = False
                        has_aryl_halide = False
                        has_benzimidazole_reactant = False

                        for reactant in reactants:
                            if checker.check_fg("Boronic acid", reactant):
                                has_boronic_acid = True
                                print(f"Found boronic acid reactant: {reactant}")
                            if checker.check_fg("Aromatic halide", reactant):
                                has_aryl_halide = True
                                print(f"Found aromatic halide reactant: {reactant}")
                            if checker.check_ring("benzimidazole", reactant):
                                has_benzimidazole_reactant = True
                                print(f"Found benzimidazole in reactant: {reactant}")

                        if (
                            has_boronic_acid
                            and has_aryl_halide
                            and has_benzimidazole_reactant
                        ):
                            late_stage_coupling = True
                            print(
                                f"Detected late-stage coupling (boronic acid + aryl halide) at depth {depth}"
                            )

            # Check for nitro reduction - early stage (higher depth)
            if depth >= 3:
                # Check for nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    nitro_reduction = True
                    print(f"Detected nitro reduction at depth {depth}")
                else:
                    # Fallback to checking reactants and products
                    has_nitro = False

                    for reactant in reactants:
                        if checker.check_fg("Nitro group", reactant):
                            has_nitro = True
                            print(f"Found nitro group in reactant: {reactant}")

                    if (
                        has_nitro
                        and checker.check_fg("Primary amine", product)
                        and not checker.check_fg("Nitro group", product)
                    ):
                        nitro_reduction = True
                        print(f"Detected nitro reduction at depth {depth}")

            # Check for SNAr reaction - early stage (higher depth)
            if depth >= 4:
                # Check for nucleophilic aromatic substitution
                if checker.check_reaction(
                    "nucl_sub_aromatic_para_nitro", rsmi
                ) or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi):
                    snar_reaction = True
                    print(
                        f"Detected nucleophilic aromatic substitution at depth {depth}"
                    )
                else:
                    # Look for F/Cl/Br displacement by NH2 (nucleophilic aromatic substitution)
                    for reactant in reactants:
                        if checker.check_fg(
                            "Aromatic halide", reactant
                        ) and checker.check_fg("Nitro group", reactant):
                            if checker.check_fg(
                                "Primary amine", product
                            ) and not checker.check_fg("Aromatic halide", product):
                                snar_reaction = True
                                print(f"Detected SNAr reaction at depth {depth}")
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Print final status
    print(f"Benzimidazole formation: {benzimidazole_formation}")
    print(f"Late-stage coupling: {late_stage_coupling}")
    print(f"Nitro reduction: {nitro_reduction}")
    print(f"SNAr reaction: {snar_reaction}")

    # Return True if all key transformations are detected
    return (
        benzimidazole_formation
        and late_stage_coupling
        and nitro_reduction
        and snar_reaction
    )
