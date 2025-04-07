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
    Detects if the synthesis route involves reduction of a nitro group followed by cyclization.
    """
    nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction

        if (
            node["type"] == "reaction" and depth <= 3
        ):  # Expanded to depth 3 to catch more potential patterns
            try:
                # Get reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for nitro reduction reaction
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                # Check for heterocycle formation reactions
                is_heterocycle_formation = checker.check_reaction(
                    "Formation of NOS Heterocycles", rsmi
                )
                is_indole_formation = checker.check_reaction("Fischer indole", rsmi)
                is_benzimidazole_formation = checker.check_reaction(
                    "benzimidazole_derivatives_aldehyde", rsmi
                ) or checker.check_reaction("benzimidazole_derivatives_carboxylic-acid/ester", rsmi)

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(reactant_mols):
                    # Check for nitro group in reactants
                    has_nitro = any(
                        [
                            checker.check_fg("Nitro group", Chem.MolToSmiles(mol))
                            for mol in reactant_mols
                        ]
                    )

                    # Check for amine in product
                    has_amine = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )

                    # Check for nitro group absence in product
                    nitro_gone = not checker.check_fg("Nitro group", product_smiles)

                    # Check for specific heterocyclic rings in reactants and product
                    heterocycles = [
                        "indole",
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                        "quinoline",
                        "isoquinoline",
                        "pyrrole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                    ]

                    reactant_rings = set()
                    for mol in reactant_mols:
                        mol_smiles = Chem.MolToSmiles(mol)
                        for ring_name in heterocycles:
                            if checker.check_ring(ring_name, mol_smiles):
                                reactant_rings.add(ring_name)

                    product_rings = set()
                    for ring_name in heterocycles:
                        if checker.check_ring(ring_name, product_smiles):
                            product_rings.add(ring_name)

                    new_heterocycle_formed = len(product_rings - reactant_rings) > 0

                    # Also check general ring count
                    reactant_ring_count = sum(
                        [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                    )
                    product_ring_count = product_mol.GetRingInfo().NumRings()
                    ring_count_increased = product_ring_count > reactant_ring_count

                    ring_formed = ring_count_increased or new_heterocycle_formed

                    print(
                        f"Nitro in reactants: {has_nitro}, Amine in product: {has_amine}, Nitro gone: {nitro_gone}"
                    )
                    print(f"Ring formed: {ring_formed}, New heterocycle: {new_heterocycle_formed}")
                    print(f"Reactant rings: {reactant_rings}, Product rings: {product_rings}")

                    # Case 1: Direct nitro reduction with cyclization
                    if is_nitro_reduction and has_nitro and has_amine and ring_formed:
                        nitro_reduction = True
                        print(f"Detected nitro reduction and cyclization at depth {depth}")

                    # Case 2: Heterocycle formation that might involve nitro reduction
                    if (
                        (
                            is_heterocycle_formation
                            or is_indole_formation
                            or is_benzimidazole_formation
                        )
                        and has_nitro
                        and nitro_gone
                    ):
                        nitro_reduction = True
                        print(
                            f"Detected heterocycle formation involving nitro group at depth {depth}"
                        )

                    # Case 3: Any reaction that reduces nitro and forms rings
                    if has_nitro and nitro_gone and ring_formed:
                        nitro_reduction = True
                        print(
                            f"Detected potential nitro reduction and cyclization in a single step at depth {depth}"
                        )

                # Case 4: Two-step process - cyclization following nitro reduction
                if depth <= 2:  # Only check for two-step process in late stages
                    # Check if this is a cyclization step
                    if product_mol and all(reactant_mols):
                        # Check for amine in reactants (potential result of previous nitro reduction)
                        has_amine_reactant = any(
                            [
                                checker.check_fg("Primary amine", Chem.MolToSmiles(mol))
                                or checker.check_fg("Secondary amine", Chem.MolToSmiles(mol))
                                or checker.check_fg("Tertiary amine", Chem.MolToSmiles(mol))
                                or checker.check_fg("Aniline", Chem.MolToSmiles(mol))
                                for mol in reactant_mols
                            ]
                        )

                        # Check for ring formation
                        reactant_ring_count = sum(
                            [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                        )
                        product_ring_count = product_mol.GetRingInfo().NumRings()
                        ring_formed = product_ring_count > reactant_ring_count

                        print(
                            f"Cyclization check at depth {depth}: Amine in reactants: {has_amine_reactant}, Ring formed: {ring_formed}"
                        )

                        if has_amine_reactant and ring_formed:
                            # Check all children for nitro reduction
                            for child in node.get("children", []):
                                if child["type"] == "reaction":
                                    try:
                                        child_rsmi = child["metadata"]["rsmi"]
                                        print(
                                            f"Checking previous step at depth {depth+1}: {child_rsmi}"
                                        )

                                        # Check for explicit nitro reduction
                                        if checker.check_reaction(
                                            "Reduction of nitro groups to amines", child_rsmi
                                        ):
                                            nitro_reduction = True
                                            print(
                                                f"Found nitro reduction in previous step: {child_rsmi}"
                                            )
                                            break

                                        # Also check for implicit nitro reduction
                                        child_reactants = child_rsmi.split(">")[0].split(".")
                                        child_product = child_rsmi.split(">")[-1]

                                        has_nitro_child = any(
                                            [
                                                checker.check_fg("Nitro group", r)
                                                for r in child_reactants
                                            ]
                                        )
                                        has_amine_child_product = (
                                            checker.check_fg("Primary amine", child_product)
                                            or checker.check_fg("Secondary amine", child_product)
                                            or checker.check_fg("Tertiary amine", child_product)
                                            or checker.check_fg("Aniline", child_product)
                                        )

                                        if has_nitro_child and has_amine_child_product:
                                            nitro_reduction = True
                                            print(
                                                f"Found implicit nitro reduction in previous step: {child_rsmi}"
                                            )
                                            break
                                    except Exception as e:
                                        print(f"Error checking child reaction: {e}")
            except Exception as e:
                print(f"Error processing node at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitro_reduction
