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
    Detects a strategy involving benzyl protection and deprotection of phenol.
    Looks for a pattern where a benzyl group is used to protect a phenol and later removed.
    """
    # Track if we found protection and deprotection
    found_benzyl_protected_phenol = False
    found_phenol_deprotection = False

    def dfs_traverse(node):
        nonlocal found_benzyl_protected_phenol, found_phenol_deprotection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzyl protected phenol pattern
                benzyl_phenol_pattern = Chem.MolFromSmarts("[c]-[O]-[CH2]-[c]")
                phenol_pattern = Chem.MolFromSmarts("[c]-[OH]")

                # Check for benzyl protected phenol in reactants or products
                for smiles in reactants + [product]:
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol and mol.HasSubstructMatch(benzyl_phenol_pattern):
                            found_benzyl_protected_phenol = True
                            print(f"Found benzyl protected phenol in: {smiles}")
                    except:
                        continue

                # Check for phenol deprotection (benzyl protected phenol → phenol)
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        any(
                            mol and mol.HasSubstructMatch(benzyl_phenol_pattern)
                            for mol in reactant_mols
                        )
                        and product_mol
                        and product_mol.HasSubstructMatch(phenol_pattern)
                        and not product_mol.HasSubstructMatch(benzyl_phenol_pattern)
                    ):
                        found_phenol_deprotection = True
                        print(f"Found phenol deprotection: {rsmi}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection are found
    return found_benzyl_protected_phenol and found_phenol_deprotection
