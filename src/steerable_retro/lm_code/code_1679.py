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
    This function detects Suzuki coupling for biaryl formation in the synthesis route.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain aryl bromide and boronic acid derivative
                has_aryl_bromide = any(
                    "[c:][Br:]" in r or "c-Br" in r or "Br[c" in r for r in reactants
                )
                has_boronic = any("B" in r and "O" in r for r in reactants)

                # Check if product contains biaryl bond that wasn't in reactants
                if has_aryl_bromide and has_boronic:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            # Look for biaryl bonds in product
                            for bond in prod_mol.GetBonds():
                                if (
                                    bond.GetBeginAtom().GetIsAromatic()
                                    and bond.GetEndAtom().GetIsAromatic()
                                    and bond.GetBeginAtom().GetAtomicNum() == 6
                                    and bond.GetEndAtom().GetAtomicNum() == 6
                                ):
                                    suzuki_coupling_found = True
                                    print("Suzuki coupling detected for biaryl formation")
                                    break
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_coupling_found
