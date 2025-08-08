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
    Detects if the route involves multiple C-N bond formations.
    """
    cn_bond_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                carboxyl_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8]")
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")

                has_amine = False
                has_carboxyl = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                        if mol.HasSubstructMatch(carboxyl_pattern):
                            has_carboxyl = True

                product_mol = Chem.MolFromSmiles(product)
                has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if has_amine and has_carboxyl and has_amide:
                    cn_bond_formations += 1

                # Check for SNAr with amine
                thioether_pattern = Chem.MolFromSmarts("[#6][#16][#6]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                cn_pattern = Chem.MolFromSmarts("[#6][#7]")

                has_thioether = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(thioether_pattern):
                            has_thioether = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                if (
                    has_thioether
                    and has_amine
                    and product_mol
                    and product_mol.HasSubstructMatch(cn_pattern)
                ):
                    cn_bond_formations += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    result = cn_bond_formations >= 2
    print(f"Multiple C-N bond formations detection: {result} (count: {cn_bond_formations})")
    return result
