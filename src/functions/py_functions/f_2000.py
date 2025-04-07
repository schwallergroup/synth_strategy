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
    Detects a synthetic strategy where a nitro group is reduced to an amine
    that participates in a late-stage lactamization to form a medium-sized ring.
    """
    # Track if we found the pattern
    found_pattern = False
    # Track if nitro group is present in early steps
    nitro_in_early_steps = False
    # Track if lactamization occurs in late steps
    lactam_formation_in_late_steps = False
    # Track if nitro reduction occurs in late steps
    nitro_reduction_in_late_steps = False

    # Track molecules with nitro groups
    nitro_containing_molecules = []
    # Track molecules with amines from nitro reduction
    amine_from_nitro_molecules = []

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, nitro_in_early_steps, lactam_formation_in_late_steps, nitro_reduction_in_late_steps
        nonlocal nitro_containing_molecules, amine_from_nitro_molecules

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for nitro group in molecules (early steps have higher depth)
            if depth >= 2:  # Early in synthesis
                if checker.check_fg("Nitro group", mol_smiles):
                    nitro_in_early_steps = True
                    nitro_containing_molecules.append(mol_smiles)
                    print(
                        f"Found nitro group in early step at depth {depth}: {mol_smiles}"
                    )

        elif node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitro reduction in any step
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction to amine at depth {depth}: {rsmi}")

                    # Check if any reactant has a nitro group that we've tracked
                    nitro_reactant = None
                    for r in reactants_smiles:
                        if checker.check_fg("Nitro group", r):
                            nitro_reactant = r
                            break

                    # Check if product has an amine
                    if nitro_reactant and (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                    ):
                        amine_from_nitro_molecules.append(product_smiles)

                        # If this is a late step, mark nitro reduction in late steps
                        if depth <= 2:
                            nitro_reduction_in_late_steps = True
                            print(
                                f"Found nitro reduction to amine in late step at depth {depth}"
                            )

                # Check for lactamization in late steps
                if depth <= 2:  # Late in synthesis
                    # Check for lactam formation reactions
                    lactam_reaction = (
                        checker.check_reaction(
                            "Intramolecular amination (heterocycle formation)", rsmi
                        )
                        or checker.check_reaction(
                            "Intramolecular transesterification/Lactone formation", rsmi
                        )
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "Ester with primary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                    )

                    # First check if any reactant has an amine that came from nitro reduction
                    amine_reactant = None
                    for r in reactants_smiles:
                        if (
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                        ) and r in amine_from_nitro_molecules:
                            amine_reactant = r
                            break

                    # Then check if product has a lactam (amide in a ring)
                    product_has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    # Check for medium-sized rings in product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        ring_info = product_mol.GetRingInfo()
                        # Check for rings of size 7-12 (medium-sized)
                        medium_rings = False
                        for ring_size in range(7, 13):
                            if any(
                                ring_info.IsAtomInRingOfSize(i, ring_size)
                                for i in range(product_mol.GetNumAtoms())
                            ):
                                medium_rings = True
                                break

                        # If product has amide and medium-sized rings, and reactant had amine from nitro
                        if (
                            amine_reactant
                            and product_has_amide
                            and medium_rings
                            and lactam_reaction
                        ):
                            lactam_formation_in_late_steps = True
                            print(
                                f"Found lactam formation in late step at depth {depth}"
                            )

                        # Also check if this is a one-pot nitro reduction and lactamization
                        elif (
                            any(
                                checker.check_fg("Nitro group", r)
                                for r in reactants_smiles
                            )
                            and product_has_amide
                            and medium_rings
                        ):
                            nitro_reduction_in_late_steps = True
                            lactam_formation_in_late_steps = True
                            print(
                                f"Found one-pot nitro reduction and lactamization in late step at depth {depth}"
                            )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Final determination based on conditions
    if nitro_in_early_steps and nitro_reduction_in_late_steps:
        found_pattern = True
        print("Detected nitro group in early steps and nitro reduction in late steps")

    if nitro_in_early_steps and lactam_formation_in_late_steps:
        found_pattern = True
        print("Detected nitro group in early steps and lactam formation in late steps")

    # If we found both nitro in early steps and either reduction or lactamization in late steps
    if found_pattern:
        print("Detected late-stage nitro-to-amine lactamization strategy")

    return found_pattern
