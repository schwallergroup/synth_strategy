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
    Detects if the synthetic route involves nucleophilic aromatic substitution
    with amine nucleophiles.
    """
    found_amine_nucleophile_snar = False

    def dfs_traverse(node):
        nonlocal found_amine_nucleophile_snar

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amine nucleophiles in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;!$(N~[!#6;!#1])]")  # Primary or secondary amine
            aromatic_pattern = Chem.MolFromSmarts("a")

            # Check if any reactant is an amine
            amine_reactant = None
            aromatic_reactant = None

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(amine_pattern):
                            amine_reactant = mol
                        if mol.HasSubstructMatch(aromatic_pattern):
                            aromatic_reactant = mol
                except:
                    continue

            # If we have both an amine and an aromatic compound as reactants
            if amine_reactant and aromatic_reactant:
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        # Check if the product has a C-N bond that wasn't in the reactants
                        c_n_bond_pattern = Chem.MolFromSmarts("a-[#7]")
                        if prod_mol.HasSubstructMatch(c_n_bond_pattern):
                            print("Found amine nucleophile in SNAr")
                            found_amine_nucleophile_snar = True
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_amine_nucleophile_snar
