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
    Detects if the synthetic route involves a cross-coupling reaction to form a biaryl bond,
    particularly looking for reactions that join two aromatic fragments where one has a halogen handle.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol is None or any(mol is None for mol in reactant_mols):
                    print("Error creating molecules from SMILES")
                    return

                # Check for halogen-containing aromatic in reactants
                halogen_aromatic = Chem.MolFromSmarts("[#6;a]-[Br,I,Cl]")
                has_halogen = any(mol.HasSubstructMatch(halogen_aromatic) for mol in reactant_mols)

                # Check for boronate in reactants (common in Suzuki couplings)
                boronate_pattern = Chem.MolFromSmarts("[#6]-B([#8])[#8]")
                has_boronate = any(mol.HasSubstructMatch(boronate_pattern) for mol in reactant_mols)

                # Check for biaryl pattern in product but not in individual reactants
                biaryl_pattern = Chem.MolFromSmarts("[#6;a]-[#6;a]")
                product_has_biaryl = product_mol.HasSubstructMatch(biaryl_pattern)

                # Check if this is a cross-coupling reaction
                if (has_halogen or has_boronate) and product_has_biaryl:
                    print(f"Found biaryl cross-coupling reaction at depth {depth}")
                    found_pattern = True
            except:
                print("Error processing reaction SMILES")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_pattern
