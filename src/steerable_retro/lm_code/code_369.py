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
    This function detects amide coupling as part of the synthetic strategy.
    """
    amide_coupling_found = False

    def dfs_traverse(node):
        nonlocal amide_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")
                amine_pattern = Chem.MolFromSmarts("[#7H2]")

                has_acid_chloride = False
                has_amine = False

                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            if r_mol.HasSubstructMatch(acid_chloride_pattern):
                                has_acid_chloride = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                    except:
                        continue

                # Check if product has amide bond
                if has_acid_chloride and has_amine:
                    try:
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol and p_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#7]-[#6](=[#8])")
                        ):
                            print("Detected amide coupling")
                            amide_coupling_found = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return amide_coupling_found
