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


def main(route):
    """
    This function detects if the synthetic route employs late-stage coupling of aromatic rings.
    """
    late_stage_coupling_found = False

    def dfs_traverse(node):
        nonlocal late_stage_coupling_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            depth = node["metadata"].get(
                "depth", 10
            )  # Default high depth if not specified

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                # Check if we're coupling aromatic rings
                if product:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Count aromatic rings in product
                        product_rings = Chem.GetSSSR(product_mol)
                        product_aromatic_rings = sum(
                            1
                            for ring in product_rings
                            if all(
                                product_mol.GetAtomWithIdx(idx).GetIsAromatic()
                                for idx in ring
                            )
                        )

                        # Count aromatic rings in each reactant
                        reactant_aromatic_rings = []
                        for reactant in reactants:
                            if reactant:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    reactant_rings = Chem.GetSSSR(reactant_mol)
                                    aromatic_count = sum(
                                        1
                                        for ring in reactant_rings
                                        if all(
                                            reactant_mol.GetAtomWithIdx(
                                                idx
                                            ).GetIsAromatic()
                                            for idx in ring
                                        )
                                    )
                                    reactant_aromatic_rings.append(aromatic_count)

                        # If product has more aromatic rings than any single reactant, rings were coupled
                        if product_aromatic_rings > max(
                            reactant_aromatic_rings, default=0
                        ):
                            print(
                                f"Late-stage aromatic coupling detected at depth {depth}"
                            )
                            late_stage_coupling_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_coupling_found
