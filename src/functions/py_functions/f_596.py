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
    This function detects a strategy involving late-stage epoxide formation
    followed by phenol O-alkylation in a linear synthetic sequence.
    """
    epoxide_molecules = {}  # Track molecules containing epoxides by their SMILES
    epoxide_formation_reactions = {}  # Track reactions that form epoxides
    phenol_alkylation_reactions = {}  # Track phenol alkylation reactions
    molecule_depths = {}  # Track all molecules and their depths

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            # Track all molecules
            mol_smiles = node["smiles"]
            molecule_depths[mol_smiles] = depth

            # Check if molecule contains an epoxide
            if checker.check_ring("oxirane", mol_smiles):
                epoxide_molecules[mol_smiles] = depth
                print(f"Molecule with epoxide found at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for epoxide formation reaction
                # In retrosynthesis, product has epoxide, reactants don't
                if checker.check_ring("oxirane", product) and not any(
                    checker.check_ring("oxirane", r) for r in reactants
                ):
                    # Check for specific epoxide formation reactions
                    if (
                        checker.check_reaction(
                            "Williamson Ether Synthesis (intra to epoxy)", rsmi
                        )
                        or checker.check_reaction("Alcohol to ether", rsmi)
                        or checker.check_reaction("Alkene to diol", rsmi)
                    ):
                        epoxide_formation_reactions[rsmi] = depth
                        print(f"Epoxide formation detected at depth {depth}: {rsmi}")

                # Check for phenol O-alkylation with epoxide
                # In retrosynthesis, reactants have epoxide and phenol, product has neither
                if any(checker.check_ring("oxirane", r) for r in reactants) and any(
                    checker.check_fg("Phenol", r) for r in reactants
                ):
                    if not checker.check_ring(
                        "oxirane", product
                    ) and not checker.check_fg("Phenol", product):
                        if (
                            checker.check_reaction("Williamson Ether Synthesis", rsmi)
                            or checker.check_reaction(
                                "Ring opening of epoxide with amine", rsmi
                            )
                            or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                        ):
                            phenol_alkylation_reactions[rsmi] = depth
                            print(
                                f"Phenol O-alkylation with epoxide detected at depth {depth}: {rsmi}"
                            )

                # Additional check for epoxide ring opening reactions
                if any(
                    checker.check_ring("oxirane", r) for r in reactants
                ) and not checker.check_ring("oxirane", product):
                    if any(
                        checker.check_fg("Phenol", r) for r in reactants
                    ) or checker.check_fg("Ether", product):
                        phenol_alkylation_reactions[rsmi] = depth
                        print(
                            f"Epoxide ring opening with phenol detected at depth {depth}: {rsmi}"
                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze results
    if epoxide_molecules:
        print(f"Found {len(epoxide_molecules)} molecules with epoxides")

        # If we have at least two epoxide-containing molecules at different depths
        if len(epoxide_molecules) >= 2:
            epoxide_depths = list(epoxide_molecules.values())
            epoxide_depths.sort()

            # Check if we have a linear sequence
            # In a linear sequence, we should have a clear path from one epoxide to another
            is_linear_sequence = False

            # Check if we have detected reactions
            if epoxide_formation_reactions and phenol_alkylation_reactions:
                min_formation_depth = min(epoxide_formation_reactions.values())
                min_alkylation_depth = min(phenol_alkylation_reactions.values())

                print(f"Minimum epoxide formation depth: {min_formation_depth}")
                print(f"Minimum phenol alkylation depth: {min_alkylation_depth}")

                # In retrosynthetic analysis, higher depth = earlier stage
                # We want epoxide formation to be later (lower depth) than alkylation (higher depth)
                if min_formation_depth < min_alkylation_depth:
                    is_linear_sequence = True
                    print(
                        "Strategy detected: Late-stage epoxide formation followed by phenol O-alkylation"
                    )
                    return True

            # If we don't have explicit reactions, infer from molecule patterns
            else:
                print(
                    "No explicit reactions detected, inferring from molecule patterns"
                )

                # Sort epoxide molecules by depth
                sorted_epoxide_mols = sorted(
                    epoxide_molecules.items(), key=lambda x: x[1]
                )

                # Check if the earliest epoxide molecule (lowest depth) has an ether but no phenol
                # This would suggest it resulted from phenol alkylation
                earliest_epoxide_smiles, earliest_depth = sorted_epoxide_mols[0]

                # Check if any later epoxide molecule has a phenol
                for later_smiles, later_depth in sorted_epoxide_mols[1:]:
                    if (
                        checker.check_fg("Phenol", later_smiles)
                        and later_depth > earliest_depth
                    ):
                        # This suggests a phenol was used with an epoxide in an earlier step
                        if checker.check_fg(
                            "Ether", earliest_epoxide_smiles
                        ) and not checker.check_fg("Phenol", earliest_epoxide_smiles):
                            print(
                                f"Inferred strategy: Late-stage epoxide formation followed by phenol O-alkylation"
                            )
                            print(
                                f"Early epoxide at depth {earliest_depth}: {earliest_epoxide_smiles}"
                            )
                            print(
                                f"Later epoxide with phenol at depth {later_depth}: {later_smiles}"
                            )
                            return True

                # Alternative inference: if we have exactly two epoxide molecules and one is at a much lower depth
                if len(epoxide_molecules) == 2 and (
                    epoxide_depths[1] - epoxide_depths[0] >= 2
                ):
                    # The deeper epoxide is likely a reagent for alkylation
                    deeper_epoxide = [
                        smiles
                        for smiles, d in epoxide_molecules.items()
                        if d == epoxide_depths[1]
                    ][0]
                    shallower_epoxide = [
                        smiles
                        for smiles, d in epoxide_molecules.items()
                        if d == epoxide_depths[0]
                    ][0]

                    # Check if the shallower epoxide has an ether linkage
                    if checker.check_fg("Ether", shallower_epoxide):
                        print(
                            f"Inferred strategy from depth difference: Late-stage epoxide formation followed by phenol O-alkylation"
                        )
                        print(
                            f"Shallow epoxide at depth {epoxide_depths[0]}: {shallower_epoxide}"
                        )
                        print(
                            f"Deep epoxide at depth {epoxide_depths[1]}: {deeper_epoxide}"
                        )
                        return True

    print("Strategy not detected")
    return False
