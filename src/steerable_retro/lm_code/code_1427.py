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
    Detects if the synthesis includes multiple C-N bond forming reactions
    (e.g., amide formation, reductive amination).
    """
    c_n_bond_formations = 0

    def dfs_traverse(node):
        nonlocal c_n_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Convert to molecules
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]
                    prod_mol = Chem.MolFromSmiles(product)

                    if not prod_mol or not reactant_mols:
                        return

                    # Check for amide formation
                    amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")
                    amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(N#*)]")

                    # Check if product has an amide that wasn't in reactants
                    if prod_mol.HasSubstructMatch(amide_pattern):
                        amide_in_reactants = False
                        for r_mol in reactant_mols:
                            if r_mol.HasSubstructMatch(amide_pattern):
                                amide_in_reactants = True
                                break

                        if not amide_in_reactants:
                            c_n_bond_formations += 1
                            print("Detected C-N bond formation: amide formation")

                    # Check for reductive amination (amine formation)
                    if prod_mol.HasSubstructMatch(amine_pattern):
                        # Check if product has a C-N bond that wasn't in reactants
                        # This is a simplified check and might need refinement
                        aldehyde_pattern = Chem.MolFromSmarts("[#6](=[#8])[#1,#6]")
                        aldehyde_in_reactants = False

                        for r_mol in reactant_mols:
                            if r_mol.HasSubstructMatch(aldehyde_pattern):
                                aldehyde_in_reactants = True
                                break

                        if aldehyde_in_reactants:
                            c_n_bond_formations += 1
                            print("Detected C-N bond formation: reductive amination")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return c_n_bond_formations >= 2
