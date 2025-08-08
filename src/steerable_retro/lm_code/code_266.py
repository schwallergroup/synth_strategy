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
    This function detects a strategy involving late-stage pyrazole formation
    from a hydrazine intermediate, with sequential nitrogen transformations.
    """
    # Track if we found key elements of the strategy
    found_pyrazole_formation = False
    found_hydrazine_intermediate = False
    found_amine_intermediate = False
    found_nitro_intermediate = False

    # Define SMARTS patterns
    pyrazole_pattern = Chem.MolFromSmarts("[n]1[n]cc[c]1")
    hydrazine_pattern = Chem.MolFromSmarts("[N][N]")
    amine_pattern = Chem.MolFromSmarts("[c][NH2]")
    nitro_pattern = Chem.MolFromSmarts("[c][N+](=[O])[O-]")

    def dfs_traverse(node, depth=0):
        nonlocal found_pyrazole_formation, found_hydrazine_intermediate
        nonlocal found_amine_intermediate, found_nitro_intermediate

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check for pyrazole formation in late stage (depth 0-1)
                if depth <= 1 and product_mol and product_mol.HasSubstructMatch(pyrazole_pattern):
                    # Check if any reactant has hydrazine
                    for reactant_smiles in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        if reactant_mol and reactant_mol.HasSubstructMatch(hydrazine_pattern):
                            found_pyrazole_formation = True
                            print(f"Found pyrazole formation at depth {depth}")

                # Check for hydrazine intermediate
                if product_mol and product_mol.HasSubstructMatch(hydrazine_pattern):
                    found_hydrazine_intermediate = True
                    print(f"Found hydrazine intermediate at depth {depth}")

                # Check for amine intermediate
                if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                    found_amine_intermediate = True
                    print(f"Found amine intermediate at depth {depth}")

                # Check for nitro intermediate
                if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                    found_nitro_intermediate = True
                    print(f"Found nitro intermediate at depth {depth}")

            except:
                print(f"Error processing SMILES at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        found_pyrazole_formation
        and found_hydrazine_intermediate
        and found_amine_intermediate
        and found_nitro_intermediate
    )

    print(f"Late-stage pyrazole formation strategy detected: {strategy_present}")
    return strategy_present
