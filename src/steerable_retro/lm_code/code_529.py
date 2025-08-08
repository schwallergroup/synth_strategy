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
    This function detects if the synthetic route involves coupling of an amino alcohol fragment.
    """
    has_amino_alcohol_coupling = False

    def dfs_traverse(node):
        nonlocal has_amino_alcohol_coupling

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if one of the reactants is an amino alcohol
                for reactant in reactants_smiles:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Pattern for amino alcohol: molecule with both NH2 and OH groups
                        amino_pattern = Chem.MolFromSmarts("[NH2]")
                        alcohol_pattern = Chem.MolFromSmarts("[OH]")

                        if mol.HasSubstructMatch(amino_pattern) and mol.HasSubstructMatch(
                            alcohol_pattern
                        ):
                            # Check if the product has both fragments combined
                            product_mol = Chem.MolFromSmiles(product_smiles)
                            if product_mol:
                                # If product has more atoms than any single reactant, it's likely a coupling
                                if product_mol.GetNumAtoms() > mol.GetNumAtoms():
                                    has_amino_alcohol_coupling = True
                                    print("Amino alcohol fragment coupling detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_amino_alcohol_coupling
