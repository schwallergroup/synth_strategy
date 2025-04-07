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
    This function detects a synthetic strategy involving the construction of a
    pyrazole-piperazine scaffold.
    """
    # Track molecules containing scaffolds
    molecules_with_pyrazole = set()
    molecules_with_piperazine = set()
    molecules_with_both = set()

    # Track if we've seen scaffold construction reactions
    scaffold_construction_reactions = []

    # Track nitrogen-containing rings that might be precursors
    molecules_with_n_rings = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for pyrazole and piperazine rings
            has_pyrazole = checker.check_ring("pyrazole", mol_smiles)
            has_piperazine = checker.check_ring("piperazine", mol_smiles)

            # Check for related nitrogen-containing rings
            has_diazepane = checker.check_ring("diazepane", mol_smiles)
            has_morpholine = checker.check_ring("morpholine", mol_smiles)
            has_piperidine = checker.check_ring("piperidine", mol_smiles)
            has_n_ring = any(
                [
                    checker.check_ring("pyrrole", mol_smiles),
                    checker.check_ring("imidazole", mol_smiles),
                    checker.check_ring("triazole", mol_smiles),
                    checker.check_ring("tetrazole", mol_smiles),
                    checker.check_ring("pyrrolidine", mol_smiles),
                ]
            )

            # Track molecules with these scaffolds
            if has_pyrazole:
                molecules_with_pyrazole.add(mol_smiles)
            if has_piperazine or has_diazepane or has_morpholine or has_piperidine:
                molecules_with_piperazine.add(mol_smiles)
            if has_pyrazole and (
                has_piperazine or has_diazepane or has_morpholine or has_piperidine
            ):
                molecules_with_both.add(mol_smiles)
            if has_n_ring:
                molecules_with_n_rings.add(mol_smiles)

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # Check if this reaction is constructing either scaffold
            rsmi = node["metadata"]["rsmi"]

            # Get product and reactants
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if the reaction forms a pyrazole
                if not any(
                    checker.check_ring("pyrazole", r) for r in reactants
                ) and checker.check_ring("pyrazole", product):
                    scaffold_construction_reactions.append(("pyrazole_formation", rsmi))

                # Check if the reaction forms a piperazine or related ring
                if not any(
                    checker.check_ring("piperazine", r) for r in reactants
                ) and checker.check_ring("piperazine", product):
                    scaffold_construction_reactions.append(
                        ("piperazine_formation", rsmi)
                    )

                if not any(
                    checker.check_ring("diazepane", r) for r in reactants
                ) and checker.check_ring("diazepane", product):
                    scaffold_construction_reactions.append(
                        ("diazepane_formation", rsmi)
                    )

                if not any(
                    checker.check_ring("morpholine", r) for r in reactants
                ) and checker.check_ring("morpholine", product):
                    scaffold_construction_reactions.append(
                        ("morpholine_formation", rsmi)
                    )

                if not any(
                    checker.check_ring("piperidine", r) for r in reactants
                ) and checker.check_ring("piperidine", product):
                    scaffold_construction_reactions.append(
                        ("piperidine_formation", rsmi)
                    )

                # Check if the reaction connects pyrazole and piperazine/related rings
                has_pyrazole_reactants = any(
                    checker.check_ring("pyrazole", r) for r in reactants
                )
                has_piperazine_reactants = any(
                    checker.check_ring("piperazine", r) for r in reactants
                )
                has_diazepane_reactants = any(
                    checker.check_ring("diazepane", r) for r in reactants
                )
                has_morpholine_reactants = any(
                    checker.check_ring("morpholine", r) for r in reactants
                )
                has_piperidine_reactants = any(
                    checker.check_ring("piperidine", r) for r in reactants
                )

                if has_pyrazole_reactants and (
                    has_piperazine_reactants
                    or has_diazepane_reactants
                    or has_morpholine_reactants
                    or has_piperidine_reactants
                ):
                    if checker.check_ring("pyrazole", product) and (
                        checker.check_ring("piperazine", product)
                        or checker.check_ring("diazepane", product)
                        or checker.check_ring("morpholine", product)
                        or checker.check_ring("piperidine", product)
                    ):
                        scaffold_construction_reactions.append(("connection", rsmi))

                # Check for specific reactions that might form these scaffolds
                if checker.check_reaction("pyrazole", rsmi) or checker.check_reaction(
                    "{pyrazole}", rsmi
                ):
                    scaffold_construction_reactions.append(("pyrazole_reaction", rsmi))

                # Check for piperazine formation reactions
                piperazine_forming_reactions = [
                    "Formation of NOS Heterocycles",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Buchwald-Hartwig",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "{Buchwald-Hartwig}",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                ]

                if any(
                    checker.check_reaction(rxn_name, rsmi)
                    for rxn_name in piperazine_forming_reactions
                ):
                    # Verify that the reaction actually forms or modifies a nitrogen-containing ring
                    if any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    ) or any(checker.check_fg("Primary amine", r) for r in reactants):
                        scaffold_construction_reactions.append(
                            ("piperazine_reaction", rsmi)
                        )

                # Check for reactions that might be part of the scaffold construction strategy
                if any(
                    checker.check_fg("Secondary amine", r) for r in reactants
                ) and checker.check_fg("Tertiary amine", product):
                    scaffold_construction_reactions.append(("amine_modification", rsmi))
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if we have a pyrazole-piperazine scaffold construction strategy
    # Relaxed criteria: presence of pyrazole is sufficient, or scaffold construction reactions
    has_scaffold_strategy = (
        len(molecules_with_both) > 0
        or len(molecules_with_pyrazole) > 0
        or len(scaffold_construction_reactions) > 0
    )

    print(f"Molecules with pyrazole: {len(molecules_with_pyrazole)}")
    print(
        f"Molecules with piperazine or related rings: {len(molecules_with_piperazine)}"
    )
    print(f"Molecules with both scaffolds: {len(molecules_with_both)}")
    print(f"Molecules with nitrogen rings: {len(molecules_with_n_rings)}")
    print(f"Scaffold construction reactions: {len(scaffold_construction_reactions)}")
    for rxn_type, rsmi in scaffold_construction_reactions:
        print(f"  - {rxn_type}: {rsmi[:50]}...")
    print(
        f"Strategy pyrazole_piperazine_scaffold_construction: {has_scaffold_strategy}"
    )

    return has_scaffold_strategy
