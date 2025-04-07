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
    Detects a late-stage convergent synthesis via N-alkylation connecting two complex fragments,
    where one fragment contains a secondary amine and the other is a chloro-substituted heterocycle.
    """
    found_pattern = False

    def is_complex_fragment(smiles):
        """Check if a molecule is a complex fragment (has multiple rings or >15 atoms)"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Count rings
            ring_info = mol.GetRingInfo()
            ring_count = ring_info.NumRings()

            # Count atoms
            atom_count = mol.GetNumAtoms()

            return ring_count > 1 or atom_count > 15
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check reactions at depths 0-2 (late-stage)
        if node["type"] == "reaction" and depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if we have at least 2 reactants (convergent synthesis)
                if len(reactants) >= 2:
                    print(
                        f"Found {len(reactants)} reactants - checking for convergent synthesis pattern"
                    )

                    # Check for any N-alkylation pattern
                    is_n_alkylation = checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    ) or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )

                    # If not directly identified, check for N-alkylation pattern by examining reactants and products
                    if not is_n_alkylation:
                        print(
                            "Not a standard N-alkylation reaction, checking for general N-alkylation pattern"
                        )
                        # Check if any reactant has an amine (primary or secondary) and another has a halide
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants
                        )
                        has_halide = any(
                            checker.check_fg("Primary halide", r)
                            or checker.check_fg("Secondary halide", r)
                            or checker.check_fg("Tertiary halide", r)
                            or checker.check_fg("Aromatic halide", r)
                            for r in reactants
                        )

                        # Check if product has a more substituted amine than reactants
                        has_more_substituted_amine = False
                        if has_amine:
                            if any(
                                checker.check_fg("Secondary amine", r)
                                for r in reactants
                            ) and checker.check_fg("Tertiary amine", product):
                                has_more_substituted_amine = True
                            elif any(
                                checker.check_fg("Primary amine", r) for r in reactants
                            ) and (
                                checker.check_fg("Secondary amine", product)
                                or checker.check_fg("Tertiary amine", product)
                            ):
                                has_more_substituted_amine = True

                        is_n_alkylation = (
                            has_amine and has_halide and has_more_substituted_amine
                        )

                        if not is_n_alkylation:
                            print("Not an N-alkylation pattern")
                            return

                    print("Confirmed N-alkylation reaction or pattern")

                    # Look for amine in one reactant and halide-heterocycle in another
                    amine_reactant = None
                    halide_heterocycle_reactant = None

                    for reactant in reactants:
                        # Check for amine
                        if checker.check_fg(
                            "Secondary amine", reactant
                        ) or checker.check_fg("Primary amine", reactant):
                            print(f"Found amine in reactant: {reactant}")
                            amine_reactant = reactant

                        # Check for halide substituent
                        has_halide = (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Aromatic halide", reactant)
                        )

                        if has_halide:
                            print(f"Found halide substituent in reactant: {reactant}")

                        # Check for any nitrogen-containing heterocycle
                        has_heterocycle = any(
                            checker.check_ring(ring, reactant)
                            for ring in [
                                "pyridine",
                                "pyrimidine",
                                "pyrazine",
                                "pyridazine",
                                "triazole",
                                "tetrazole",
                                "imidazole",
                                "oxazole",
                                "thiazole",
                                "pyrrole",
                                "indole",
                                "quinoline",
                                "isoquinoline",
                                "purine",
                                "carbazole",
                                "acridine",
                                "benzimidazole",
                                "benzoxazole",
                                "benzothiazole",
                                "piperidine",
                                "piperazine",
                                "morpholine",
                                "thiomorpholine",
                            ]
                        )

                        if has_heterocycle:
                            print(f"Found heterocycle in reactant: {reactant}")

                        if has_halide and has_heterocycle:
                            print(
                                f"Found halide-substituted heterocycle in reactant: {reactant}"
                            )
                            halide_heterocycle_reactant = reactant

                    # Check if we found both required components
                    if amine_reactant and halide_heterocycle_reactant:
                        # Check if both fragments are complex
                        amine_complex = is_complex_fragment(amine_reactant)
                        heterocycle_complex = is_complex_fragment(
                            halide_heterocycle_reactant
                        )

                        print(f"Amine fragment complexity: {amine_complex}")
                        print(f"Heterocycle fragment complexity: {heterocycle_complex}")

                        # Check if product has a more substituted amine (result of N-alkylation)
                        has_n_alkylation_product = False
                        if checker.check_fg("Tertiary amine", product):
                            print(f"Found tertiary amine in product: {product}")
                            has_n_alkylation_product = True
                        elif checker.check_fg("Secondary amine", product) and any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        ):
                            print(
                                f"Found secondary amine in product (from primary amine): {product}"
                            )
                            has_n_alkylation_product = True

                        # For convergent synthesis, both fragments should be complex
                        if (
                            has_n_alkylation_product
                            and amine_complex
                            and heterocycle_complex
                        ):
                            found_pattern = True
                            print("Found late-stage convergent N-alkylation pattern")
                        else:
                            if not has_n_alkylation_product:
                                print(
                                    "Product does not contain expected N-alkylation result"
                                )
                            else:
                                print(
                                    "Not a convergent synthesis - one or both fragments are not complex"
                                )
                    else:
                        if not amine_reactant:
                            print("Missing amine reactant")
                        if not halide_heterocycle_reactant:
                            print("Missing halide-heterocycle reactant")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_pattern
