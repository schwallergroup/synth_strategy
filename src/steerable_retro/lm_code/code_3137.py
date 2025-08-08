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
    Detects synthesis routes that use organometallic reagents (particularly Grignard)
    to form ketones, especially from Weinreb amides or other activated carboxylic derivatives.
    """
    found_organometallic_ketone_formation = False

    def dfs_traverse(node):
        nonlocal found_organometallic_ketone_formation

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for organometallic reagent (Grignard or organolithium)
            grignard_pattern = Chem.MolFromSmarts("[Mg][#6]")
            organolithium_pattern = Chem.MolFromSmarts("[Li][#6]")

            # Check for ketone in product
            ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")

            # Check reactants for organometallic reagent
            has_organometallic = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(grignard_pattern) or mol.HasSubstructMatch(
                            organolithium_pattern
                        ):
                            has_organometallic = True
                            break
                        # Also check for text indicators of Grignard in SMILES
                        if "[Mg]" in reactant or "MgCl" in reactant or "MgBr" in reactant:
                            has_organometallic = True
                            break
                except:
                    # If SMILES parsing fails, check for text indicators
                    if "[Mg]" in reactant or "MgCl" in reactant or "MgBr" in reactant:
                        has_organometallic = True
                        break

            # Check product for ketone
            try:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(ketone_pattern) and has_organometallic:
                    found_organometallic_ketone_formation = True
                    print("Found organometallic ketone formation")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_organometallic_ketone_formation
