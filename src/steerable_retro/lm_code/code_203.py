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
    Detects synthesis strategies that use bromo groups as coupling handles
    to form C-C bonds between aromatic fragments.
    """
    has_bromo_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_bromo_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for bromo in reactants
            bromo_reactants = []
            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("c[Br]")):
                    bromo_reactants.append(reactant)

            # If we have bromo reactants, check if product has new C-C bond
            if bromo_reactants:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    # This is a simplification - in a real implementation we would need to
                    # analyze the reaction more carefully to confirm C-C bond formation
                    if (
                        not product_mol.HasSubstructMatch(Chem.MolFromSmarts("c[Br]"))
                        or product_mol.GetNumAtoms()
                        > sum(Chem.MolFromSmiles(r).GetNumAtoms() for r in reactants_smiles) - 5
                    ):
                        has_bromo_coupling = True
                        print(f"Detected bromo-to-CC coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_bromo_coupling
