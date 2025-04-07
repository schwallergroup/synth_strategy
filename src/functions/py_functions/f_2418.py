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
    Detects a convergent synthesis with multiple ring formations and fragment couplings.
    """
    ring_formations = 0
    fragment_couplings = 0
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal ring_formations, fragment_couplings, late_stage_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Count reactant fragments
                num_reactants = len(reactants_smiles)

                # Check for ring formation
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if (
                        all(mol is not None for mol in reactant_mols)
                        and product_mol is not None
                    ):
                        reactant_ring_count = sum(
                            Chem.GetSSSR(mol) for mol in reactant_mols
                        )
                        product_ring_count = Chem.GetSSSR(product_mol)

                        if product_ring_count > reactant_ring_count:
                            ring_formations += 1
                            print(f"Ring formation detected at depth {depth}")

                        # Check for fragment coupling (multiple reactants to single product)
                        if num_reactants > 1:
                            fragment_couplings += 1
                            print(f"Fragment coupling detected at depth {depth}")

                            # Check if it's late-stage (depth < 3)
                            if depth < 3:
                                late_stage_coupling = True
                                print(
                                    f"Late-stage fragment coupling detected at depth {depth}"
                                )
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if we have multiple ring formations and fragment couplings
    return ring_formations >= 1 and fragment_couplings >= 1 and late_stage_coupling
