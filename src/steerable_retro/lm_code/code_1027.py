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
    Detects lactam formation as a key ring-forming strategy.
    """
    lactam_formation_found = False

    def dfs_traverse(node):
        nonlocal lactam_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amine in reactants
            amine_in_reactants = False
            for r in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r)
                if r_mol:
                    amine_pattern = Chem.MolFromSmarts("[#7;!$(N~[!#6]);!$(N~[#6]=[#8])]")
                    if r_mol.HasSubstructMatch(amine_pattern):
                        amine_in_reactants = True
                        break

            # Check for lactam in product
            if amine_in_reactants:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    lactam_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])")
                    if product_mol.HasSubstructMatch(lactam_pattern):
                        # Check if this is a ring formation
                        for r in reactants_smiles:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol and r_mol.HasSubstructMatch(lactam_pattern):
                                # Lactam already present in reactant
                                break
                        else:
                            # Lactam not found in any reactant
                            lactam_formation_found = True
                            print(f"Found lactam formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return lactam_formation_found
