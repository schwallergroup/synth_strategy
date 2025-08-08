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
    This function detects the use of Boc protection throughout synthesis
    with removal in the final steps.
    """
    boc_present = False
    boc_removed = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_present, boc_removed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Boc group SMARTS pattern
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)")

                # Check for Boc in reactants
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                        boc_present = True

                # Check for Boc removal (present in reactant but not in product)
                reactant_has_boc = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                        reactant_has_boc = True

                product_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_boc
                    and product_mol
                    and not product_mol.HasSubstructMatch(boc_pattern)
                    and depth <= 1
                ):
                    boc_removed = True
                    print(f"Detected Boc removal at depth {depth}: {rsmi}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return boc_present and boc_removed
