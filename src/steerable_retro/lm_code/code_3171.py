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
    Detects if the synthesis route uses a sequence of nucleophilic substitutions,
    particularly looking for halogen replacement by N or S nucleophiles.
    """
    substitution_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for nucleophilic substitution
                c_x_pattern = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")  # C-Cl, C-Br, C-I
                c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")  # C-N
                c_s_pattern = Chem.MolFromSmarts("[#6]-[#16]")  # C-S

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check each reactant for halogen pattern
                    for reactant in reactants:
                        if not reactant:
                            continue

                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        # If reactant has C-X bond
                        if reactant_mol.HasSubstructMatch(c_x_pattern):
                            # Check if product has new C-N or C-S bonds
                            if product_mol.HasSubstructMatch(c_n_pattern):
                                print(f"Nucleophilic substitution: C-X to C-N at depth {depth}")
                                substitution_reactions.append((depth, "C-N"))

                            if product_mol.HasSubstructMatch(c_s_pattern):
                                print(f"Nucleophilic substitution: C-X to C-S at depth {depth}")
                                substitution_reactions.append((depth, "C-S"))
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have at least 2 nucleophilic substitution reactions
    return len(substitution_reactions) >= 2
