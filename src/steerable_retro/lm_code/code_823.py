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
    This function detects if the synthetic route involves N-alkylation of a heterocycle (e.g., indole).
    """
    has_n_alkylation = False

    def dfs_traverse(node):
        nonlocal has_n_alkylation

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                # Check for N-methylated heterocycle in product
                n_methyl_heterocycle = Chem.MolFromSmarts(
                    "[#6]1:[#6]:[n]([CH3]):[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:12"
                )

                if product_mol.HasSubstructMatch(n_methyl_heterocycle):
                    # Check if one reactant is an NH heterocycle
                    nh_heterocycle = Chem.MolFromSmarts(
                        "[#6]1:[#6]:[nH]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:12"
                    )

                    # Check if another reactant is a methyl source (e.g., MeI)
                    methyl_source = Chem.MolFromSmarts("C[I,Br,Cl,O]")

                    reactant_has_nh = False
                    reactant_has_methyl = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        if reactant_mol.HasSubstructMatch(nh_heterocycle):
                            reactant_has_nh = True
                        if reactant_mol.HasSubstructMatch(methyl_source):
                            reactant_has_methyl = True

                    if reactant_has_nh and reactant_has_methyl:
                        has_n_alkylation = True
                        print(f"Detected N-alkylation of heterocycle: {rsmi}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_n_alkylation
