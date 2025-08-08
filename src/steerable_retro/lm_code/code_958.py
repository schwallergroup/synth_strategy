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
    Detects nucleophilic aromatic substitution of chlorine with an amine.
    """
    nas_found = False

    def dfs_traverse(node):
        nonlocal nas_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nucleophilic aromatic substitution
                try:
                    # Look for aryl chloride in reactants
                    aryl_chloride_pattern = Chem.MolFromSmarts("c-[Cl]")
                    amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                    has_aryl_chloride = False
                    has_amine = False

                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if not react_mol:
                            continue

                        if react_mol.HasSubstructMatch(aryl_chloride_pattern):
                            has_aryl_chloride = True

                        if react_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                    # Check for C-N bond formation in product
                    if has_aryl_chloride and has_amine:
                        aryl_amine_pattern = Chem.MolFromSmarts("c-[#7;!$(N-[!#6;!#1])]")
                        prod_mol = Chem.MolFromSmiles(product)

                        if prod_mol and prod_mol.HasSubstructMatch(aryl_amine_pattern):
                            # Check if chlorine count decreased
                            cl_count_reactants = sum(
                                [reactant.count("Cl") for reactant in reactants]
                            )
                            cl_count_product = product.count("Cl")

                            if cl_count_product < cl_count_reactants:
                                nas_found = True
                                print("Nucleophilic aromatic substitution detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nas_found
