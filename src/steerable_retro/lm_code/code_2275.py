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


def main(route):
    """
    This function detects if the synthesis follows a linear pattern
    (each step builds on a single previous product).

    A linear synthesis is defined as a route where:
    1. Each reaction node has at most one reaction child
    2. The main structural backbone is preserved throughout the synthesis
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # Check if this is a reaction node
        if node["type"] == "reaction" and "children" in node:
            # Get product and reactants
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants_smiles = reactants_part.split(".")
                product_mol = Chem.MolFromSmiles(product_part)

                # Count reaction children (non-leaf mol nodes that lead to reactions)
                reaction_children = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        if any(
                            grandchild["type"] == "reaction"
                            for grandchild in child.get("children", [])
                        ):
                            reaction_children += 1

                # If more than one reaction child, it's not linear
                if reaction_children > 1:
                    print(f"Non-linear synthesis detected (branching at depth {depth}): {rsmi}")
                    is_linear = False
                    return

                # If we have multiple reactants, check if one is clearly the main contributor
                if len(reactants_smiles) > 1 and product_mol:
                    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                    reactant_mols = [mol for mol in reactant_mols if mol]  # Filter out None values

                    # Skip if we couldn't parse molecules
                    if not reactant_mols:
                        return

                    # Find the reactant with the largest common substructure with the product
                    max_mcs_size = 0
                    for reactant_mol in reactant_mols:
                        if (
                            reactant_mol.GetNumAtoms() < 3
                        ):  # Skip very small molecules (likely reagents)
                            continue

                        mcs = rdFMCS.FindMCS(
                            [reactant_mol, product_mol],
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            bondCompare=rdFMCS.BondCompare.CompareOrder,
                            completeRingsOnly=True,
                            timeout=1,
                        )

                        if mcs.numAtoms > max_mcs_size:
                            max_mcs_size = mcs.numAtoms

                    # If no significant common substructure with any reactant, might be convergent
                    product_size = product_mol.GetNumAtoms()
                    if max_mcs_size < product_size * 0.5:  # Less than 50% similarity
                        print(
                            f"Potentially non-linear synthesis detected (low structural similarity): {rsmi}"
                        )
                        # We don't set is_linear=False here as this is just a warning
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return is_linear
