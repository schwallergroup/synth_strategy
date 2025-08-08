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
    Detects if the synthesis involves a protection-deprotection strategy,
    particularly focusing on Cbz protection of amines.
    """
    has_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Cbz protection pattern
                # Amine reactant + Cbz-Cl → Cbz-protected amine
                amine_pattern = Chem.MolFromSmarts("[NH2,NH]")
                cbz_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1[CH2][O][C](=[O])[N]")
                cbz_reagent_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#6][#6][#6][#6]1[CH2][O][C](=[O])[Cl,OH]"
                )

                has_amine_reactant = False
                has_cbz_reagent = False
                has_protected_product = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(amine_pattern):
                            has_amine_reactant = True
                        if mol and mol.HasSubstructMatch(cbz_reagent_pattern):
                            has_cbz_reagent = True
                    except:
                        continue

                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    has_protected_product = prod_mol and prod_mol.HasSubstructMatch(cbz_pattern)
                except:
                    has_protected_product = False

                if has_amine_reactant and (has_cbz_reagent or has_protected_product):
                    has_protection = True
                    print(f"Detected Cbz protection strategy at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_protection
