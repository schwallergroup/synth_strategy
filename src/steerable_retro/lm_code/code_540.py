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
    This function detects if the synthesis route involves phenol O-alkylation
    with an alpha-carbonyl electrophile.
    """
    phenol_alkylation_found = False

    def dfs_traverse(node):
        nonlocal phenol_alkylation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for phenol in reactants
            phenol_pattern = Chem.MolFromSmarts("c[OX2H1]")

            # Check for alpha-carbonyl electrophile in reactants
            alpha_carbonyl_pattern = Chem.MolFromSmarts("[Br,Cl,I,F][CX4][CX3](=[OX1])")

            # Check for phenoxy-carbonyl in product
            phenoxy_carbonyl_pattern = Chem.MolFromSmarts("c[OX2][CX4][CX3](=[OX1])")

            phenol_found = False
            alpha_carbonyl_found = False

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                if reactant_mol.HasSubstructMatch(phenol_pattern):
                    phenol_found = True

                if reactant_mol.HasSubstructMatch(alpha_carbonyl_pattern):
                    alpha_carbonyl_found = True

            product_mol = Chem.MolFromSmiles(product_smiles)
            product_has_phenoxy_carbonyl = product_mol and product_mol.HasSubstructMatch(
                phenoxy_carbonyl_pattern
            )

            if phenol_found and alpha_carbonyl_found and product_has_phenoxy_carbonyl:
                print("Phenol O-alkylation detected in reaction:", rsmi)
                phenol_alkylation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return phenol_alkylation_found
