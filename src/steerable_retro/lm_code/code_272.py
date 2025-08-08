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
    Detects if the synthesis involves benzylation of an amine with a benzyl halide.
    """
    found_benzylation = False

    def dfs_traverse(node):
        nonlocal found_benzylation

        if node["type"] == "reaction":
            # Extract reactants and product
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product and reactants:
                    # Check for benzyl halide and amine in reactants
                    benzyl_halide_pattern = Chem.MolFromSmarts("[Br,Cl,I,F][CH2]c1ccccc1")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    # Check for benzylated amine in product
                    benzylated_amine_pattern = Chem.MolFromSmarts("[NH]([CH2]c1ccccc1)")

                    has_benzyl_halide = any(
                        r and r.HasSubstructMatch(benzyl_halide_pattern) for r in reactants
                    )
                    has_amine = any(r and r.HasSubstructMatch(amine_pattern) for r in reactants)
                    has_benzylated_amine = product and product.HasSubstructMatch(
                        benzylated_amine_pattern
                    )

                    if has_benzyl_halide and has_amine and has_benzylated_amine:
                        found_benzylation = True
                        print("Found benzylation of amine with benzyl halide")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_benzylation
