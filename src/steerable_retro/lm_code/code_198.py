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
    This function detects if the synthesis route involves a late-stage N-alkylation
    with a fluorobenzyl group.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            # Check if this is an N-alkylation reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for fluorobenzyl pattern in reactants
                fluorobenzyl_pattern = Chem.MolFromSmarts("[#6]-[#6]1[#6][#6][#6]([#9])[#6][#6]1")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(fluorobenzyl_pattern):
                            # Check if the other reactant has an NH group
                            for other_reactant in reactants:
                                if other_reactant != reactant:
                                    other_mol = Chem.MolFromSmiles(other_reactant)
                                    nh_pattern = Chem.MolFromSmarts("[N;H]")
                                    if other_mol and other_mol.HasSubstructMatch(nh_pattern):
                                        # Check if product has the fluorobenzyl attached to N
                                        prod_mol = Chem.MolFromSmiles(product)
                                        if prod_mol:
                                            result = True
                                            print(
                                                f"Found late-stage N-alkylation with fluorobenzyl at depth {depth}"
                                            )
                    except:
                        continue

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return result
