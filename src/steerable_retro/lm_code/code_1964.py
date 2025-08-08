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
    This function detects if the synthesis utilizes halogen chemistry extensively
    (multiple halogenated intermediates or halogenation reactions).
    """
    halogen_pattern = Chem.MolFromSmarts("[#6]-[#9,#17,#35,#53]")
    halogen_reactions_count = 0
    halogenated_intermediates_count = 0

    def dfs_traverse(node):
        nonlocal halogen_reactions_count, halogenated_intermediates_count

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(halogen_pattern):
                halogenated_intermediates_count += 1

        elif node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
            ]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check if this reaction introduces or removes a halogen
                product_halogens = len(product_mol.GetSubstructMatches(halogen_pattern))
                reactant_halogens = sum(
                    len(r.GetSubstructMatches(halogen_pattern)) for r in reactant_mols if r
                )

                if product_halogens != reactant_halogens:
                    halogen_reactions_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Consider it a halogen chemistry strategy if we have multiple halogenated intermediates
    # and at least one halogenation/dehalogenation reaction
    result = halogenated_intermediates_count >= 3 and halogen_reactions_count >= 1
    print(
        f"Halogen chemistry strategy detected: {result} (intermediates: {halogenated_intermediates_count}, reactions: {halogen_reactions_count})"
    )
    return result
