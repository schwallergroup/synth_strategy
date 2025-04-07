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
    Detects a late-stage N-alkylation of a piperazine with a benzofuran-containing fragment.
    """
    has_piperazine = False
    has_benzofuran = False
    has_alkylation = False

    def is_benzofuran_or_derivative(smiles):
        """Check if molecule contains benzofuran core structure"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            # Check for benzofuran ring directly
            if checker.check_ring("benzofuran", smiles):
                return True

            # Check for benzofuran-like structure using substructure match
            benzofuran_pattern = Chem.MolFromSmarts("c1cccc2occ(*)c12")
            if mol and benzofuran_pattern and mol.HasSubstructMatch(benzofuran_pattern):
                print(f"Found benzofuran-like structure in: {smiles}")
                return True

            return False
        except:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal has_piperazine, has_benzofuran, has_alkylation

        # Update node depth for tracking
        node["depth"] = depth

        if node["type"] == "reaction" and depth <= 2:  # Check first three steps (late-stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for piperazine in reactants
                piperazine_reactant = None
                has_leaving_group = False
                benzofuran_reactant = None

                # First check if any reactant has piperazine
                for reactant in reactants:
                    if checker.check_ring("piperazine", reactant):
                        has_piperazine = True
                        piperazine_reactant = reactant
                        print(f"Found piperazine in reactant: {reactant}")

                # Then check if any reactant has benzofuran or derivative
                for reactant in reactants:
                    if is_benzofuran_or_derivative(reactant):
                        has_benzofuran = True
                        benzofuran_reactant = reactant
                        print(f"Found benzofuran or derivative in reactant: {reactant}")

                # Check for leaving groups in any reactant
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary halide", reactant)
                        or checker.check_fg("Secondary halide", reactant)
                        or checker.check_fg("Tertiary halide", reactant)
                        or checker.check_fg("Mesylate", reactant)
                        or checker.check_fg("Tosylate", reactant)
                        or checker.check_fg("Triflate", reactant)
                    ):
                        has_leaving_group = True
                        print(f"Found leaving group in reactant: {reactant}")

                # Check if this is an N-alkylation reaction
                if piperazine_reactant and benzofuran_reactant and has_leaving_group:
                    # Check for various alkylation reaction types
                    if (
                        checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    ):

                        product_mol = Chem.MolFromSmiles(product)
                        if (
                            product_mol
                            and checker.check_ring("piperazine", product)
                            and is_benzofuran_or_derivative(product)
                        ):
                            has_alkylation = True
                            print(f"Found late-stage piperazine alkylation at depth {depth}")

                # If we haven't found the reaction yet, try a more general approach
                if not has_alkylation and piperazine_reactant and has_leaving_group:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        product_mol
                        and checker.check_ring("piperazine", product)
                        and is_benzofuran_or_derivative(product)
                    ):
                        # If product contains both piperazine and benzofuran, and we had the right reactants,
                        # it's likely an alkylation even if not explicitly matched by reaction type
                        has_alkylation = True
                        print(
                            f"Found likely late-stage piperazine alkylation at depth {depth} (general check)"
                        )

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = has_piperazine and has_benzofuran and has_alkylation
    print(f"Late-stage piperazine alkylation strategy detected: {result}")
    print(f"  - Piperazine present: {has_piperazine}")
    print(f"  - Benzofuran present: {has_benzofuran}")
    print(f"  - Alkylation reaction: {has_alkylation}")

    return result
