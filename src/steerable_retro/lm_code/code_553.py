#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    Detects a synthetic strategy involving late-stage N-alkylation with halogenated intermediates,
    connecting complex heterocyclic fragments.
    """
    n_alkylation_found = False
    ketone_to_alcohol_found = False

    # List of heterocyclic rings to check
    heterocyclic_rings = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "furan",
        "thiophene",
    ]

    def has_heterocycle(smiles):
        """Check if a molecule contains heterocyclic rings"""
        for ring in heterocyclic_rings:
            if checker.check_ring(ring, smiles):
                print(f"Found heterocycle {ring} in {smiles}")
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_found, ketone_to_alcohol_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for N-alkylation at late stage (depth 0-1)
                if depth <= 1:
                    # Check for N-alkylation reaction
                    is_n_alkylation = checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    ) or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )

                    # Fallback to checking functional groups if reaction check fails
                    if not is_n_alkylation:
                        # Check for halides and amines in reactants
                        has_halide = any(
                            checker.check_fg("Primary halide", r)
                            or checker.check_fg("Secondary halide", r)
                            or checker.check_fg("Tertiary halide", r)
                            for r in reactants_smiles
                        )

                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants_smiles
                        )

                        # Check if product has a tertiary amine (result of N-alkylation)
                        has_tertiary_amine_product = checker.check_fg(
                            "Tertiary amine", product_smiles
                        )

                        is_n_alkylation = has_halide and has_amine and has_tertiary_amine_product

                    # Verify that heterocyclic fragments are being connected
                    if is_n_alkylation:
                        # Check if reactants contain heterocycles
                        heterocycles_in_reactants = sum(
                            has_heterocycle(r) for r in reactants_smiles
                        )

                        # Check if product contains heterocycles
                        heterocycles_in_product = has_heterocycle(product_smiles)

                        # N-alkylation should connect heterocyclic fragments
                        if heterocycles_in_reactants >= 1 and heterocycles_in_product:
                            n_alkylation_found = True
                            print(
                                f"Found late-stage N-alkylation connecting heterocycles at depth {depth}: {rsmi}"
                            )

                # Check for ketone to alcohol reduction at early-mid stage (depth 2-3)
                if 2 <= depth <= 3:
                    # Check for reduction reaction
                    is_reduction = checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    )

                    # Fallback to checking functional groups if reaction check fails
                    if not is_reduction:
                        # Check for ketone in reactants and alcohol in product
                        has_ketone = any(checker.check_fg("Ketone", r) for r in reactants_smiles)

                        has_alcohol = (
                            checker.check_fg("Primary alcohol", product_smiles)
                            or checker.check_fg("Secondary alcohol", product_smiles)
                            or checker.check_fg("Tertiary alcohol", product_smiles)
                        )

                        is_reduction = has_ketone and has_alcohol

                    # Verify that the reduction is relevant to a heterocyclic system
                    if is_reduction:
                        # Check if reactants or product contain heterocycles
                        if any(has_heterocycle(r) for r in reactants_smiles) or has_heterocycle(
                            product_smiles
                        ):
                            ketone_to_alcohol_found = True
                            print(
                                f"Found ketone to alcohol reduction related to heterocycles at depth {depth}: {rsmi}"
                            )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if both key transformations are found
    print(
        f"N-alkylation found: {n_alkylation_found}, Ketone to alcohol found: {ketone_to_alcohol_found}"
    )
    return n_alkylation_found and ketone_to_alcohol_found
