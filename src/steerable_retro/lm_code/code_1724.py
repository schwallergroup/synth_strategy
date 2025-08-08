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
    Detects a linear synthetic strategy where fragments are assembled
    sequentially through amide bond formations.
    """
    amide_coupling_steps = []
    linear_assembly = True

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_steps, linear_assembly

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is potentially an amide coupling
                if len(reactants_smiles) >= 2:  # At least two reactants
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                    product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                    if product and all(r is not None for r in reactants):
                        # Check for amide formation
                        carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                        amine_pattern = Chem.MolFromSmarts("[NH2]")
                        amide_pattern = Chem.MolFromSmarts("[C](=O)[NH]")

                        if (
                            any(r.HasSubstructMatch(carboxylic_acid_pattern) for r in reactants)
                            and any(r.HasSubstructMatch(amine_pattern) for r in reactants)
                            and product.HasSubstructMatch(amide_pattern)
                        ):
                            amide_coupling_steps.append(depth)
                            print(f"Found amide coupling at depth {depth}")

                            # Check if more than 2 reactants - would indicate convergent rather than linear
                            if len(reactants) > 2:
                                linear_assembly = False
                                print(
                                    f"Found convergent step with {len(reactants)} reactants at depth {depth}"
                                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort amide coupling steps by depth
    amide_coupling_steps.sort()

    # Check if we have multiple amide couplings in a linear fashion
    has_multiple_amide_couplings = len(amide_coupling_steps) >= 2

    print(f"Amide coupling steps: {amide_coupling_steps}")
    print(f"Linear assembly: {linear_assembly}")

    return has_multiple_amide_couplings and linear_assembly
