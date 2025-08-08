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
    Detects if the synthetic route involves silyl protection of alcohols.
    Looks for TBDMS or TBS protection patterns.
    """
    protection_count = 0

    def dfs_traverse(node):
        nonlocal protection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # SMARTS for silyl ether (TBDMS/TBS)
                silyl_pattern = Chem.MolFromSmarts("[C][Si]([C])([C])[O][C]")

                # Check if product contains silyl ether but reactants don't
                if product_mol and silyl_pattern:
                    product_matches = product_mol.GetSubstructMatches(silyl_pattern)
                    if product_matches:
                        # Check if this is a protection reaction (reactants have OH)
                        alcohol_pattern = Chem.MolFromSmarts("[O][C]")
                        for r_mol in reactant_mols:
                            if (
                                r_mol
                                and alcohol_pattern
                                and r_mol.HasSubstructMatch(alcohol_pattern)
                            ):
                                # Check if reactant doesn't already have silyl group
                                if not r_mol.HasSubstructMatch(silyl_pattern):
                                    protection_count += 1
                                    print(f"Found silyl protection reaction: {rsmi}")
                                    break
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least one silyl protection is found
    return protection_count >= 1
