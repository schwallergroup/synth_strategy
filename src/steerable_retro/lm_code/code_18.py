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
    Detects if the synthetic route involves formation of a pyrazolone ring
    through reaction of a β-ketoester with hydrazine.
    """
    pyrazolone_pattern = Chem.MolFromSmarts("[#6]1=[#7][#7][#6](=[#8])[#6]1")
    hydrazine_pattern = Chem.MolFromSmarts("[NH2][NH2]")
    ketoester_pattern = Chem.MolFromSmarts("[#6][#6](=[#8])[#6][#6](=[#8])[#8][#6]")

    found_pyrazolone_formation = False

    def dfs_traverse(node):
        nonlocal found_pyrazolone_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if hydrazine is a reactant
                hydrazine_found = False
                ketoester_found = False

                for reactant in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(hydrazine_pattern):
                            hydrazine_found = True
                        if mol and mol.HasSubstructMatch(ketoester_pattern):
                            ketoester_found = True
                    except:
                        continue

                # Check if product has pyrazolone
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(pyrazolone_pattern):
                        if hydrazine_found and ketoester_found:
                            print("Found pyrazolone formation from ketoester and hydrazine")
                            found_pyrazolone_formation = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_pyrazolone_formation
