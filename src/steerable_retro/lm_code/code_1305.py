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
    This function detects multiple functional group interconversions,
    such as alcohol to bromide, alcohol to carboxylic acid, etc.
    """
    interconversions = []

    def dfs_traverse(node):
        nonlocal interconversions

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and reactant_mols:
                    # Alcohol to bromide
                    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                    bromide_pattern = Chem.MolFromSmarts("[#6]-[#35]")

                    if any(
                        r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols if r
                    ) and product_mol.HasSubstructMatch(bromide_pattern):
                        interconversions.append("alcohol_to_bromide")
                        print(f"Found alcohol to bromide conversion: {rsmi}")

                    # Benzyl alcohol to carboxylic acid
                    benzyl_alcohol_pattern = Chem.MolFromSmarts("[c]-[#6]-[#8;H1]")
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[c]-[#6](=[#8])-[#8;H1]")

                    if any(
                        r.HasSubstructMatch(benzyl_alcohol_pattern) for r in reactant_mols if r
                    ) and product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                        interconversions.append("benzyl_alcohol_to_carboxylic_acid")
                        print(f"Found benzyl alcohol to carboxylic acid conversion: {rsmi}")

                    # Aryl halide to phenol
                    aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")
                    phenol_pattern = Chem.MolFromSmarts("[c]-[#8;H1]")

                    if any(
                        r.HasSubstructMatch(aryl_halide_pattern) for r in reactant_mols if r
                    ) and product_mol.HasSubstructMatch(phenol_pattern):
                        interconversions.append("aryl_halide_to_phenol")
                        print(f"Found aryl halide to phenol conversion: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    unique_interconversions = set(interconversions)
    print(f"Unique functional group interconversions: {unique_interconversions}")

    # Return True if we have at least 2 different types of interconversions
    return len(unique_interconversions) >= 2
