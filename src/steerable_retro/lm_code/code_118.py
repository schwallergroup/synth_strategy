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
    Detects if stereocenters are preserved throughout the synthesis.
    """
    has_stereocenter = False
    stereocenters_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal has_stereocenter, stereocenters_preserved

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(r for r in reactants):
                # Check for stereocenters in product
                product_chiral_centers = Chem.FindMolChiralCenters(product, includeUnassigned=False)

                if product_chiral_centers:
                    has_stereocenter = True

                    # Check if all stereocenters in reactants are preserved
                    reactant_chiral_centers = []
                    for r in reactants:
                        reactant_chiral_centers.extend(
                            Chem.FindMolChiralCenters(r, includeUnassigned=False)
                        )

                    # Simple check - if number of stereocenters changes, they're not preserved
                    if len(product_chiral_centers) != len(reactant_chiral_centers):
                        stereocenters_preserved = False

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if has_stereocenter and stereocenters_preserved:
        print("Detected stereocenter preservation strategy")
        return True
    else:
        print("Stereocenter preservation strategy not detected")
        return False
