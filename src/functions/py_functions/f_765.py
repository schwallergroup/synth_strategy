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
    This function detects if the synthetic route involves sequential nitrogen functionalization steps.
    """
    # Track molecules that have undergone N-functionalization
    n_functionalized_molecules = {}
    # Track sequential steps
    sequential_steps = False

    def dfs_traverse(node, depth=0):
        nonlocal sequential_steps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-functionalization reactions
                is_n_functionalization = False

                # Check for common N-functionalization reactions
                if (
                    checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        rsmi,
                    )
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )
                ):

                    is_n_functionalization = True
                    print(f"Found N-functionalization reaction: {rsmi}")

                # If not detected by reaction type, check for amine reactants and N-containing products
                if not is_n_functionalization:
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )

                    if has_amine and (
                        checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Sulfonamide", product)
                        or checker.check_fg("Urea", product)
                        or checker.check_fg("Thiourea", product)
                        or checker.check_fg("Carbamic ester", product)
                    ):
                        is_n_functionalization = True
                        print(f"Found N-functionalization by FG analysis: {rsmi}")

                # If this is an N-functionalization, track the product
                if is_n_functionalization:
                    # Store the product as having undergone N-functionalization
                    # Convert to canonical SMILES for consistent comparison
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        canonical_product = Chem.MolToSmiles(product_mol)
                        n_functionalized_molecules[canonical_product] = depth
                    else:
                        n_functionalized_molecules[product] = depth

                    # Check if any reactant has previously undergone N-functionalization
                    for reactant in reactants:
                        # Convert reactant to canonical SMILES
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            canonical_reactant = Chem.MolToSmiles(reactant_mol)
                            # Check if this reactant is in our tracked molecules
                            if canonical_reactant in n_functionalized_molecules:
                                print(
                                    f"Sequential N-functionalization detected! Previous at depth {n_functionalized_molecules[canonical_reactant]}, current at depth {depth}"
                                )
                                sequential_steps = True
                            else:
                                # Also check by substructure match for atom-mapped molecules
                                for (
                                    processed_smiles,
                                    processed_depth,
                                ) in n_functionalized_molecules.items():
                                    processed_mol = Chem.MolFromSmiles(processed_smiles)
                                    if processed_mol and (
                                        Chem.MolFromSmiles(reactant).HasSubstructMatch(
                                            processed_mol
                                        )
                                        or processed_mol.HasSubstructMatch(
                                            Chem.MolFromSmiles(reactant)
                                        )
                                    ):
                                        print(
                                            f"Sequential N-functionalization detected via substructure! Previous at depth {processed_depth}, current at depth {depth}"
                                        )
                                        sequential_steps = True
                                        break
                        else:
                            # Fallback to direct string comparison if SMILES parsing fails
                            if reactant in n_functionalized_molecules:
                                print(
                                    f"Sequential N-functionalization detected (direct match)! Previous at depth {n_functionalized_molecules[reactant]}, current at depth {depth}"
                                )
                                sequential_steps = True

        # Process children (going backward in synthesis)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Count the number of N-functionalized molecules
    n_functionalization_steps = len(n_functionalized_molecules)
    print(f"Total N-functionalization steps: {n_functionalization_steps}")
    print(f"Sequential steps detected: {sequential_steps}")

    # Return True if we have at least 2 N-functionalization steps and they're sequential
    return n_functionalization_steps >= 2 and sequential_steps
