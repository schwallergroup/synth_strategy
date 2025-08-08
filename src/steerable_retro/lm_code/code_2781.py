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
    This function detects if the synthetic route employs a late-stage nitro reduction strategy,
    where a nitro group is carried through multiple steps before being reduced to an amine
    in the final or near-final step.
    """
    # Track if we found a nitro reduction
    nitro_reduction_found = False
    # Track the depth at which nitro reduction occurs
    nitro_reduction_depth = None
    # Track if nitro group exists in earlier steps
    nitro_in_earlier_steps = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, nitro_reduction_depth, nitro_in_earlier_steps

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and any(reactants):
                # Check for nitro groups in reactants
                nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
                reactants_with_nitro = any(
                    mol.HasSubstructMatch(nitro_pattern) for mol in reactants if mol
                )

                # Check for amine groups in product
                amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
                product_has_amine = product.HasSubstructMatch(amine_pattern) if product else False

                # Check if this reaction is a nitro reduction
                if (
                    reactants_with_nitro and product_has_amine and depth <= 1
                ):  # Depth <= 1 means late stage
                    print(f"Found nitro reduction at depth {depth}")
                    nitro_reduction_found = True
                    nitro_reduction_depth = depth

                # Check if nitro exists in earlier steps
                if reactants_with_nitro and depth > 1:
                    print(f"Found nitro group at earlier depth {depth}")
                    nitro_in_earlier_steps = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found a nitro reduction at late stage
    # and the nitro group was present in earlier steps
    strategy_present = nitro_reduction_found and nitro_in_earlier_steps

    print(f"Late-stage nitro reduction strategy detected: {strategy_present}")
    if strategy_present:
        print(f"Nitro reduction occurred at depth {nitro_reduction_depth}")

    return strategy_present
