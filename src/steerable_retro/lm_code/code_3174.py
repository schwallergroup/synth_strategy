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
    This function detects the combined strategy of cyclopropanation with Weinreb amide
    intermediates and late-stage amide coupling.
    """
    # Initialize flags for each strategy component
    has_cyclopropane = False
    has_weinreb_amide = False
    has_late_amide_coupling = False
    has_amide_in_target = False

    # Track depth for late-stage determination
    max_depth = 0
    amide_coupling_depth = None

    # Check if the target molecule (root node) contains an amide group
    if route["type"] == "mol":
        target_smiles = route["smiles"]
        if (
            checker.check_fg("Primary amide", target_smiles)
            or checker.check_fg("Secondary amide", target_smiles)
            or checker.check_fg("Tertiary amide", target_smiles)
        ):
            print(f"Target molecule contains amide group: {target_smiles}")
            has_amide_in_target = True

    def dfs_traverse(node, depth=0):
        nonlocal has_cyclopropane, has_weinreb_amide, has_late_amide_coupling, max_depth, amide_coupling_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        # Process molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for cyclopropane ring
            if checker.check_ring("cyclopropane", mol_smiles):
                print(f"Found cyclopropane ring in molecule: {mol_smiles}")
                has_cyclopropane = True

            # Check for Weinreb amide functional group
            # Weinreb amides have N-methoxy-N-methyl amide structure
            if "CON(C)C(=O)" in mol_smiles or "C(=O)N(C)OC" in mol_smiles:
                print(f"Found potential Weinreb amide in molecule: {mol_smiles}")
                # Verify it's actually a Weinreb amide by checking the structure
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "N" and atom.GetDegree() == 3:
                            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                            if "C" in neighbors and "O" in neighbors:
                                carbon_neighbors = [
                                    n for n in atom.GetNeighbors() if n.GetSymbol() == "C"
                                ]
                                for c in carbon_neighbors:
                                    if any(
                                        n.GetSymbol() == "O" and n.GetDegree() == 1
                                        for n in c.GetNeighbors()
                                    ):
                                        has_weinreb_amide = True
                                        print(f"Confirmed Weinreb amide in molecule: {mol_smiles}")
                                        break

        # Process reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for cyclopropanation reaction
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check if product has cyclopropane but reactants don't
            if checker.check_ring("cyclopropane", product) and not any(
                checker.check_ring("cyclopropane", r) for r in reactants
            ):
                print(f"Found cyclopropanation reaction: {rxn_smiles}")
                has_cyclopropane = True

            # Check for amide coupling reaction - expanded to catch more cases
            is_amide_coupling = False

            # Check named reactions
            if (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rxn_smiles,
                )
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    rxn_smiles,
                )
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rxn_smiles
                )
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rxn_smiles
                )
                or checker.check_reaction("Acyl chloride with secondary amine to amide", rxn_smiles)
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rxn_smiles)
                or checker.check_reaction("Ester with primary amine to amide", rxn_smiles)
                or checker.check_reaction("Ester with secondary amine to amide", rxn_smiles)
                or checker.check_reaction("Schotten-Baumann to ester", rxn_smiles)
                or checker.check_reaction("Schotten-Baumann_amide", rxn_smiles)
                or checker.check_reaction("Acylation of primary amines", rxn_smiles)
                or checker.check_reaction("Acylation of secondary amines", rxn_smiles)
            ):
                is_amide_coupling = True

            # Check for amide formation by functional group appearance
            elif (
                (
                    checker.check_fg("Primary amide", product)
                    and not any(checker.check_fg("Primary amide", r) for r in reactants)
                )
                or (
                    checker.check_fg("Secondary amide", product)
                    and not any(checker.check_fg("Secondary amide", r) for r in reactants)
                )
                or (
                    checker.check_fg("Tertiary amide", product)
                    and not any(checker.check_fg("Tertiary amide", r) for r in reactants)
                )
            ):
                # Verify it's an amide formation by checking for amine and carboxylic acid/acyl halide reactants
                has_amine = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    or checker.check_fg("Aniline", r)
                    for r in reactants
                )
                has_acyl = any(
                    checker.check_fg("Carboxylic acid", r)
                    or checker.check_fg("Acyl halide", r)
                    or checker.check_fg("Ester", r)
                    or checker.check_fg("Anhydride", r)
                    for r in reactants
                )
                if has_amine and has_acyl:
                    is_amide_coupling = True

            if is_amide_coupling:
                print(f"Found amide coupling reaction at depth {depth}: {rxn_smiles}")
                # If this is the first amide coupling or it's at a lower depth (later stage)
                if amide_coupling_depth is None or depth < amide_coupling_depth:
                    amide_coupling_depth = depth

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Determine if amide coupling is late-stage (lower depth)
    if amide_coupling_depth is not None:
        # Consider it late-stage if it's in the first 75% of the synthesis depth
        if amide_coupling_depth <= max_depth * 0.75:
            has_late_amide_coupling = True
            print(
                f"Amide coupling at depth {amide_coupling_depth} is considered late-stage (max depth: {max_depth})"
            )
    # If we didn't find an amide coupling reaction but the target has an amide, consider it late-stage
    elif has_amide_in_target:
        print("No explicit amide coupling reaction found, but target molecule contains amide group")
        has_late_amide_coupling = True

    # Check if all three strategy components are present
    combined_strategy = has_cyclopropane and has_weinreb_amide and has_late_amide_coupling

    if combined_strategy:
        print(
            "Detected combined strategy: cyclopropanation with Weinreb amide intermediates and late-stage amide coupling"
        )
    else:
        print(
            f"Strategy components: cyclopropane={has_cyclopropane}, Weinreb amide={has_weinreb_amide}, late amide coupling={has_late_amide_coupling}"
        )

    return combined_strategy
