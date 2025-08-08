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
    Detects the use of a Grignard reaction to introduce an alkyl chain to an aromatic ring,
    specifically looking for replacement of aryl-Cl with aryl-alkyl.
    """
    grignard_found = False

    def dfs_traverse(node):
        nonlocal grignard_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant contains "MgCl" or similar Grignard reagent pattern
            grignard_reactant = False
            for reactant in reactants:
                if "[Mg]" in reactant and "Cl" in reactant:
                    grignard_reactant = True
                    break

            if grignard_reactant:
                # Check for aryl-Cl to aryl-alkyl transformation
                aryl_chloride_pattern = Chem.MolFromSmarts("c-Cl")

                # Find the non-Grignard reactant
                for reactant in reactants:
                    if "[Mg]" not in reactant:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(aryl_chloride_pattern):
                            product_mol = Chem.MolFromSmiles(product)

                            # Check if product has alkyl chain where chloride was
                            if product_mol:
                                alkyl_chain_pattern = Chem.MolFromSmarts("c-C-C-C")
                                if product_mol.HasSubstructMatch(alkyl_chain_pattern):
                                    grignard_found = True
                                    print("Found Grignard reaction introducing alkyl chain")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return grignard_found
