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
    This function detects if the synthetic route involves nitro group reduction to amine.
    """
    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")
    amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
    nitro_to_amine_found = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                reactant_has_nitro = False
                for r_smiles in reactants_smiles:
                    try:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(nitro_pattern):
                            reactant_has_nitro = True
                            break
                    except:
                        continue

                # Check for amine in product
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    product_has_amine = product_mol and product_mol.HasSubstructMatch(amine_pattern)

                    if reactant_has_nitro and product_has_amine:
                        print("Nitro to amine reduction detected")
                        nitro_to_amine_found = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_to_amine_found
