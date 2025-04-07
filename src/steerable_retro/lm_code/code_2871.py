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
    Detects if the synthetic route follows a linear fragment assembly strategy
    where fragments are added sequentially rather than in a convergent manner.
    """
    reaction_count = 0
    linear_assembly_reactions = 0

    # Expanded list of coupling reactions used in fragment assembly
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with sulfonic esters",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl OTf",
        "Sonogashira acetylene_aryl OTf",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl halide",
        "Heck terminal vinyl",
        "Heck_terminal_vinyl",
        "Heck_non-terminal_vinyl",
        "Oxidative Heck reaction",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_aryl OTf",
        "Stille reaction_vinyl OTf",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Ullmann condensation",
        "Williamson Ether Synthesis",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "Chan-Lam alcohol",
        "Chan-Lam amine",
        "Chan-Lam etherification",
        "decarboxylative_coupling",
    ]

    # Track the depth of each linear assembly reaction
    linear_assembly_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_assembly_reactions

        if node["type"] == "reaction":
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a typical fragment assembly reaction
            is_linear_assembly = False

            # Condition 1: Has exactly 2 reactants (most common for fragment assembly)
            if len(reactants) == 2:
                # Condition 2: Check if it's a known coupling reaction
                for reaction_type in coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_linear_assembly = True
                        print(f"Found coupling reaction: {reaction_type}")
                        break

                # Condition 3: If not a known coupling, check for functional groups involved in coupling
                if not is_linear_assembly:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Check if reactants have complementary functional groups for coupling
                    if product_mol and all(reactant_mols):
                        # Check for boronic acid/ester and halide/triflate (Suzuki-like)
                        has_boronic = any(
                            checker.check_fg("Boronic acid", r)
                            or checker.check_fg("Boronic ester", r)
                            for r in reactants
                        )
                        has_halide = any(
                            checker.check_fg("Aromatic halide", r)
                            or checker.check_fg("Triflate", r)
                            for r in reactants
                        )

                        # Check for alkyne and halide (Sonogashira-like)
                        has_alkyne = any(checker.check_fg("Alkyne", r) for r in reactants)

                        # Check for amine and halide (Buchwald-Hartwig-like)
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants
                        )

                        if (
                            (has_boronic and has_halide)
                            or (has_alkyne and has_halide)
                            or (has_amine and has_halide)
                        ):
                            is_linear_assembly = True
                            print(f"Found complementary functional groups for coupling")

                # Condition 4: Check if product is more complex than reactants
                if not is_linear_assembly:
                    product_mol = (
                        Chem.MolFromSmiles(product)
                        if not "product_mol" in locals()
                        else product_mol
                    )
                    reactant_mols = (
                        [Chem.MolFromSmiles(r) for r in reactants]
                        if not "reactant_mols" in locals()
                        else reactant_mols
                    )

                    if product_mol and all(reactant_mols):
                        # Compare by atom count
                        product_atoms = product_mol.GetNumAtoms()
                        max_reactant_atoms = max(mol.GetNumAtoms() for mol in reactant_mols)

                        # Compare by ring count
                        product_rings = product_mol.GetRingInfo().NumRings()
                        max_reactant_rings = max(
                            mol.GetRingInfo().NumRings() for mol in reactant_mols
                        )

                        # Product should be significantly larger than largest reactant
                        # OR have more rings (indicating fragment joining)
                        if (
                            product_atoms > max_reactant_atoms * 1.2
                            or product_rings > max_reactant_rings
                        ):
                            is_linear_assembly = True
                            print(
                                f"Found fragment joining: product size {product_atoms} vs largest reactant {max_reactant_atoms}"
                            )
                            print(
                                f"Ring count: product {product_rings} vs largest reactant {max_reactant_rings}"
                            )

            if is_linear_assembly:
                linear_assembly_reactions += 1
                linear_assembly_depths.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if most reactions follow linear assembly pattern
    if reaction_count > 0:
        linear_ratio = linear_assembly_reactions / reaction_count

        print(
            f"Linear assembly ratio: {linear_ratio:.2f} ({linear_assembly_reactions}/{reaction_count} reactions are linear assembly)"
        )

        # Adjust threshold based on route length
        threshold = 0.4 if reaction_count <= 5 else 0.5

        # Check if linear assembly reactions are sequential (depths are consecutive)
        sequential_pattern = False
        if linear_assembly_depths:
            # Sort depths (remember lower depth = later stage in synthesis)
            sorted_depths = sorted(linear_assembly_depths)
            # Check if depths form a consecutive sequence
            if len(sorted_depths) >= 2:
                gaps = [
                    sorted_depths[i + 1] - sorted_depths[i] for i in range(len(sorted_depths) - 1)
                ]
                # If most gaps are small (1, 2, or 3), it's likely sequential
                sequential_pattern = sum(1 for g in gaps if g <= 3) >= len(gaps) * 0.5
                print(f"Reaction depths: {sorted_depths}, Sequential pattern: {sequential_pattern}")

        # If ratio exceeds threshold or we have sequential pattern with at least 2 linear reactions
        # or we simply have at least 2 linear assembly reactions in a short route
        if (
            linear_ratio >= threshold
            or (sequential_pattern and linear_assembly_reactions >= 2)
            or linear_assembly_reactions >= 2
        ):
            print("Linear fragment assembly strategy detected")
            return True

    print("No linear fragment assembly strategy detected")
    return False
