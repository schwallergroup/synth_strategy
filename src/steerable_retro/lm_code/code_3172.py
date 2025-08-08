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
    Detects if the synthesis route employs halogen chemistry as key intermediates,
    looking for halogenation and subsequent replacement of halogens.
    """
    halogenation_steps = []
    halogen_replacement_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for halogen chemistry
                c_x_pattern = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")  # C-Cl, C-Br, C-I

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    # Count halogen bonds in product and reactants
                    if product_mol:
                        product_c_x_count = len(product_mol.GetSubstructMatches(c_x_pattern))
                        reactant_c_x_count = sum(
                            len(mol.GetSubstructMatches(c_x_pattern))
                            for mol in reactant_mols
                            if mol
                        )

                        # If product has more halogen bonds than reactants combined (halogenation)
                        if product_c_x_count > reactant_c_x_count:
                            print(f"Halogenation detected at depth {depth}")
                            halogenation_steps.append(depth)

                        # If product has fewer halogen bonds than reactants combined (halogen replacement)
                        if product_c_x_count < reactant_c_x_count:
                            print(f"Halogen replacement detected at depth {depth}")
                            halogen_replacement_steps.append(depth)
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both halogenation and subsequent replacement
    return len(halogenation_steps) > 0 and len(halogen_replacement_steps) > 0
