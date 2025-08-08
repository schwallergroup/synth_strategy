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
    This function detects a strategy involving sequential heteroatom alkylations
    (O-methylation followed by N-methylation) with nitro group introduction early in the synthesis.
    """
    # Initialize tracking variables
    has_n_methylation = False
    has_o_methylation = False
    has_nitro_introduction = False
    n_methylation_depth = None
    o_methylation_depth = None
    nitro_introduction_depth = None

    # Define SMARTS patterns
    methylamine_pattern = Chem.MolFromSmarts("[CH3][NH2]")
    diazomethane_pattern = Chem.MolFromSmarts("[N-]=[N+]=[CH2]")
    nitro_pattern = Chem.MolFromSmarts("[#8]=[N+]([#8-])-[#6]")

    def dfs_traverse(node, depth=0):
        nonlocal has_n_methylation, has_o_methylation, has_nitro_introduction
        nonlocal n_methylation_depth, o_methylation_depth, nitro_introduction_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for N-methylation
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(methylamine_pattern):
                    has_n_methylation = True
                    n_methylation_depth = depth
                    print(f"Found N-methylation at depth {depth}")

            # Check for O-methylation
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(diazomethane_pattern):
                    has_o_methylation = True
                    o_methylation_depth = depth
                    print(f"Found O-methylation at depth {depth}")

            # Check for nitro group introduction
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                # Check if any reactant doesn't have the nitro group
                nitro_introduced = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and not reactant_mol.HasSubstructMatch(nitro_pattern):
                        nitro_introduced = True

                if nitro_introduced:
                    has_nitro_introduction = True
                    nitro_introduction_depth = depth
                    print(f"Found nitro group introduction at depth {depth}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present
    # 1. All three transformations must be present
    # 2. N-methylation should be at a lower depth than O-methylation (later in synthesis)
    # 3. Nitro introduction should be at a higher depth (earlier in synthesis)
    strategy_present = (
        has_n_methylation
        and has_o_methylation
        and has_nitro_introduction
        and n_methylation_depth < o_methylation_depth
        and o_methylation_depth < nitro_introduction_depth
    )

    print(f"Sequential heteroatom alkylation strategy detected: {strategy_present}")
    return strategy_present
