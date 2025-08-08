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
    This function detects if the synthesis includes an aromatic nucleophilic substitution
    step (converting Cl to NH2 on an aromatic ring).
    """
    snAr_found = False

    def dfs_traverse(node):
        nonlocal snAr_found

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Pattern for aromatic chloride
            aryl_chloride_pattern = Chem.MolFromSmarts("[c][Cl]")
            # Pattern for aromatic amine
            aryl_amine_pattern = Chem.MolFromSmarts("[c][NH2]")

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                # Check for SNAr
                if (
                    any(
                        r and r.HasSubstructMatch(aryl_chloride_pattern) for r in reactant_mols if r
                    )
                    and product_mol
                    and product_mol.HasSubstructMatch(aryl_amine_pattern)
                ):
                    snAr_found = True
                    print("Detected aromatic nucleophilic substitution step")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return snAr_found
