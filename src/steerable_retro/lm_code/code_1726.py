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
    Detects a synthetic strategy involving lactone ring opening followed by
    alcohol protection with silyl groups.
    """
    # Track if we found the key transformations
    found_lactone_opening = False
    found_silyl_protection = False

    def dfs_traverse(node):
        nonlocal found_lactone_opening, found_silyl_protection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not product or not all(reactants):
                return

            # Check for lactone opening (lactone in reactant, diol in product)
            lactone_pattern = Chem.MolFromSmarts("[#6]-[#8]-C(=O)-[#6]")
            diol_pattern = Chem.MolFromSmarts("[#6]-[OH].[#6]-[OH]")

            for reactant in reactants:
                if reactant.HasSubstructMatch(lactone_pattern):
                    # Check if product has diol pattern or two separate OH groups
                    if (
                        product
                        and sum(
                            1 for match in product.GetSubstructMatches(Chem.MolFromSmarts("[OH]"))
                        )
                        >= 2
                    ):
                        print("Found lactone opening to diol")
                        found_lactone_opening = True

            # Check for silyl protection (alcohol in reactant, silyl ether in product)
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH]")
            silyl_ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[Si]")

            for reactant in reactants:
                if (
                    reactant.HasSubstructMatch(alcohol_pattern)
                    and product
                    and product.HasSubstructMatch(silyl_ether_pattern)
                ):
                    print("Found silyl protection of alcohol")
                    found_silyl_protection = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both key transformations were found
    return found_lactone_opening and found_silyl_protection
