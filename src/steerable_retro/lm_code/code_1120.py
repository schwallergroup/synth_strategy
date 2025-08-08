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
    This function detects if the synthetic route involves coupling with a piperazine moiety.
    """
    piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([C])CC1")
    piperazine_coupling_found = False

    def dfs_traverse(node):
        nonlocal piperazine_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if piperazine is in reactants
                piperazine_in_reactants = False
                for r_smiles in reactants_smiles:
                    try:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(piperazine_pattern):
                            piperazine_in_reactants = True
                            break
                    except:
                        continue

                # Check if piperazine is in product and was incorporated
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if (
                        product_mol
                        and product_mol.HasSubstructMatch(piperazine_pattern)
                        and piperazine_in_reactants
                    ):
                        # Check if this is a coupling reaction (more than one reactant)
                        if len(reactants_smiles) > 1:
                            print("Piperazine coupling detected")
                            piperazine_coupling_found = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return piperazine_coupling_found
