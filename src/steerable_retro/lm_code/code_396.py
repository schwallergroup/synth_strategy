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
    This function detects early-stage amide formation from an amine and acid chloride.
    Early-stage means the amide formation occurs at a high depth in the synthetic tree.
    """
    amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected

        if node["type"] == "reaction" and depth >= 2:  # Early step (depth >= 2)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                # Check for acid chloride in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#6]")

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    has_amine = False
                    has_acid_chloride = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                            if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                                has_acid_chloride = True

                    if (
                        has_amine
                        and has_acid_chloride
                        and product_mol
                        and product_mol.HasSubstructMatch(amide_pattern)
                    ):
                        print(f"Early-stage amide formation detected: {rsmi}")
                        amide_formation_detected = True
                except:
                    print("Error processing SMILES in amide formation detection")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return amide_formation_detected
