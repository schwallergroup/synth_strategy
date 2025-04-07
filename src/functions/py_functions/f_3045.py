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
    This function detects if the synthetic route employs a late-stage SNAr reaction
    to introduce a cyclic amine onto a chloropyrimidine or similar scaffold.
    """
    result = False

    # List of cyclic amines to check
    cyclic_amine_rings = [
        "pyrrolidine",
        "piperidine",
        "azetidine",
        "azepane",
        "morpholine",
        "piperazine",
        "thiomorpholine",
    ]

    # List of heteroaromatic rings
    heteroaromatic_rings = [
        "pyrimidine",
        "pyridine",
        "pyrazine",
        "pyridazine",
        "quinoline",
        "isoquinoline",
        "purine",
        "triazole",
        "tetrazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal result

        # We're looking for a reaction at depth 0 or 1 (final or penultimate step)
        if node["type"] == "reaction" and depth <= 1:
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an SNAr or related N-arylation reaction
                is_snar = (
                    checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        rsmi,
                    )
                )

                # If not a known SNAr reaction, check manually for the pattern
                if not is_snar:
                    # Look for a cyclic amine and a halogenated heteroaromatic in reactants
                    # and check if they combine in the product
                    cyclic_amine_reactant = None
                    halogen_hetero_reactant = None

                    # First, identify the cyclic amine reactant
                    for reactant in reactants:
                        for ring in cyclic_amine_rings:
                            if checker.check_ring(ring, reactant):
                                cyclic_amine_reactant = reactant
                                print(
                                    f"Found cyclic amine ({ring}) in reactant: {reactant}"
                                )
                                break
                        if cyclic_amine_reactant:
                            break

                    # Then, identify a halogenated heteroaromatic reactant
                    for reactant in reactants:
                        if reactant == cyclic_amine_reactant:
                            continue

                        hetero_ring_found = False
                        for ring in heteroaromatic_rings:
                            if checker.check_ring(ring, reactant):
                                hetero_ring_found = True
                                print(
                                    f"Found heteroaromatic ring ({ring}) in reactant: {reactant}"
                                )
                                break

                        # Check for halogens
                        has_halogen = False
                        if (
                            "Cl" in reactant
                            or "Br" in reactant
                            or "I" in reactant
                            or "F" in reactant
                        ):
                            has_halogen = True
                        elif checker.check_fg("Aromatic halide", reactant):
                            has_halogen = True

                        if hetero_ring_found and has_halogen:
                            halogen_hetero_reactant = reactant
                            print(
                                f"Found halogenated heteroaromatic in reactant: {reactant}"
                            )
                            break

                    # If we found both components, check if they combine in the product
                    if cyclic_amine_reactant and halogen_hetero_reactant:
                        # Check if the product contains both structures
                        product_has_cyclic_amine = False
                        product_has_heteroaromatic = False

                        for ring in cyclic_amine_rings:
                            if checker.check_ring(ring, product):
                                product_has_cyclic_amine = True
                                print(
                                    f"Product contains cyclic amine ({ring}): {product}"
                                )
                                break

                        for ring in heteroaromatic_rings:
                            if checker.check_ring(ring, product):
                                product_has_heteroaromatic = True
                                print(
                                    f"Product contains heteroaromatic ring ({ring}): {product}"
                                )
                                break

                        # Check if the product has more atoms than either reactant alone
                        # This suggests they've combined
                        try:
                            product_mol = Chem.MolFromSmiles(product)
                            cyclic_amine_mol = Chem.MolFromSmiles(cyclic_amine_reactant)
                            hetero_mol = Chem.MolFromSmiles(halogen_hetero_reactant)

                            if (
                                product_mol
                                and cyclic_amine_mol
                                and hetero_mol
                                and product_mol.GetNumAtoms()
                                > cyclic_amine_mol.GetNumAtoms()
                                and product_mol.GetNumAtoms() > hetero_mol.GetNumAtoms()
                            ):

                                # Check if the halogen count decreased
                                halogen_count_reactant = (
                                    halogen_hetero_reactant.count("Cl")
                                    + halogen_hetero_reactant.count("Br")
                                    + halogen_hetero_reactant.count("I")
                                    + halogen_hetero_reactant.count("F")
                                )

                                halogen_count_product = (
                                    product.count("Cl")
                                    + product.count("Br")
                                    + product.count("I")
                                    + product.count("F")
                                )

                                if (
                                    halogen_count_product < halogen_count_reactant
                                    and product_has_cyclic_amine
                                    and product_has_heteroaromatic
                                ):
                                    print(
                                        f"Detected SNAr pattern: cyclic amine replaced halogen on heteroaromatic"
                                    )
                                    is_snar = True
                        except Exception as e:
                            print(f"Error analyzing molecules: {e}")

                if is_snar:
                    print(f"Detected SNAr reaction at depth {depth}")

                    # Check for cyclic amine in reactants
                    cyclic_amine_found = False
                    for reactant in reactants:
                        for ring in cyclic_amine_rings:
                            if checker.check_ring(ring, reactant):
                                cyclic_amine_found = True
                                print(
                                    f"Found cyclic amine ({ring}) in reactant: {reactant}"
                                )
                                break
                        if cyclic_amine_found:
                            break

                    # Check for heteroaromatic ring in reactants
                    hetero_ring_found = False
                    for reactant in reactants:
                        for ring in heteroaromatic_rings:
                            if checker.check_ring(ring, reactant):
                                hetero_ring_found = True
                                print(
                                    f"Found heteroaromatic ring ({ring}) in reactant: {reactant}"
                                )
                                break
                        if hetero_ring_found:
                            break

                    # Check if the product contains both structures
                    if cyclic_amine_found and hetero_ring_found:
                        product_has_cyclic_amine = False
                        product_has_heteroaromatic = False

                        for ring in cyclic_amine_rings:
                            if checker.check_ring(ring, product):
                                product_has_cyclic_amine = True
                                print(
                                    f"Product contains cyclic amine ({ring}): {product}"
                                )
                                break

                        for ring in heteroaromatic_rings:
                            if checker.check_ring(ring, product):
                                product_has_heteroaromatic = True
                                print(
                                    f"Product contains heteroaromatic ring ({ring}): {product}"
                                )
                                break

                        if product_has_cyclic_amine and product_has_heteroaromatic:
                            print(f"Confirmed SNAr with cyclic amine at depth {depth}")
                            result = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return result
