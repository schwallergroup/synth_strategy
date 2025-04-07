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
    This function detects aromatic bromination, typically using NBS or similar reagents.
    """
    from rdkit import Chem

    has_aromatic_bromination = False

    def dfs_traverse(node):
        nonlocal has_aromatic_bromination

        if has_aromatic_bromination:
            return  # Early return if we already found what we're looking for

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Direct check for aromatic bromination reaction types
            if checker.check_reaction("Aromatic bromination", rsmi):
                print(f"Detected aromatic bromination reaction: {rsmi}")
                has_aromatic_bromination = True
                return

            # Check for other specific bromination reactions that might be relevant
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has aromatic halide specifically bromine
            if checker.check_fg("Aromatic halide", product):
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check if the aromatic halide is specifically bromine
                    aryl_bromide_pattern = Chem.MolFromSmarts("c[Br]")
                    if product_mol.HasSubstructMatch(aryl_bromide_pattern):
                        # Check if any reactant has a brominating agent
                        has_brominating_agent = False
                        for reactant in reactants:
                            # Check for common brominating agents
                            if (
                                (
                                    "Br" in reactant
                                    and not checker.check_fg(
                                        "Aromatic halide", reactant
                                    )
                                )
                                or "NBS" in reactant.upper()
                                or "N-BROMOSUCCINIMIDE" in reactant.upper()
                            ):
                                has_brominating_agent = True
                                print(f"Found brominating agent: {reactant}")
                                break

                        if has_brominating_agent:
                            # Count aromatic bromines in reactants and product
                            product_aryl_br_count = len(
                                product_mol.GetSubstructMatches(aryl_bromide_pattern)
                            )

                            reactant_aryl_br_count = 0
                            for r in reactants:
                                # Check if reactant has aromatic structure
                                is_aromatic = False
                                for ring in [
                                    "benzene",
                                    "naphthalene",
                                    "anthracene",
                                    "pyridine",
                                    "furan",
                                    "thiophene",
                                    "pyrrole",
                                ]:
                                    if checker.check_ring(ring, r):
                                        is_aromatic = True
                                        break

                                if is_aromatic or checker.check_fg(
                                    "Aromatic halide", r
                                ):
                                    r_mol = Chem.MolFromSmiles(r)
                                    if r_mol:
                                        reactant_aryl_br_count += len(
                                            r_mol.GetSubstructMatches(
                                                aryl_bromide_pattern
                                            )
                                        )

                            if product_aryl_br_count > reactant_aryl_br_count:
                                print(
                                    f"Detected aromatic bromination by bromine count: {rsmi}"
                                )
                                has_aromatic_bromination = True
                                return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_aromatic_bromination
