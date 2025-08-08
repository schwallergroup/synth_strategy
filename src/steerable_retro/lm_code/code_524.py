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
    This function detects a synthetic strategy involving amide bond formation
    between a carboxylic acid and an amine.
    """
    amide_coupling_found = False

    def dfs_traverse(node):
        nonlocal amide_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("c-[NH2]")
                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("C(=O)N")

                has_acid = False
                has_amine = False

                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            if r_mol.HasSubstructMatch(carboxylic_acid_pattern):
                                has_acid = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                    except:
                        continue

                try:
                    p_mol = Chem.MolFromSmiles(product)
                    if has_acid and has_amine and p_mol and p_mol.HasSubstructMatch(amide_pattern):
                        print("Found amide coupling")
                        amide_coupling_found = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_coupling_found
