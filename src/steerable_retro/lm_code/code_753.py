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
    This function detects if the synthetic route involves nitration of an aromatic ring followed by reduction to amine.
    """
    nitration_step = False
    nitro_reduction_step = False

    def dfs_traverse(node):
        nonlocal nitration_step, nitro_reduction_step

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitration (adding a nitro group to an aromatic ring)
            nitro_pattern = Chem.MolFromSmarts("c[N+](=[O])[O-]")

            # Check if product has nitro group but reactants don't
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                has_nitro_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                        has_nitro_in_reactants = True
                        break

                if not has_nitro_in_reactants:
                    nitration_step = True
                    print(f"Nitration detected in reaction: {rsmi}")

            # Check for nitro reduction to amine
            amine_pattern = Chem.MolFromSmarts("c[NH2]")

            # Check if reactants have nitro group and product has amine
            has_nitro_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                    has_nitro_in_reactants = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            if (
                has_nitro_in_reactants
                and product_mol
                and product_mol.HasSubstructMatch(amine_pattern)
            ):
                nitro_reduction_step = True
                print(f"Nitro reduction detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if both nitration and reduction steps are found
    result = nitration_step and nitro_reduction_step
    print(f"Nitration-reduction strategy: {result}")
    return result
