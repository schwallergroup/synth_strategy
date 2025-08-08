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
    This function detects a cascade of protecting group manipulations:
    ethyl ester → free acid → p-nitrobenzyl ester
    """
    # Track the protecting group transformations
    ethyl_ester_deprotection = False
    pnb_ester_protection = False

    def dfs_traverse(node):
        nonlocal ethyl_ester_deprotection, pnb_ester_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Create RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Patterns for detecting protecting groups
                ethyl_ester_pattern = Chem.MolFromSmarts("C(=O)OCC")
                carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
                pnb_ester_pattern = Chem.MolFromSmarts("C(=O)OCc1ccc([N+](=O)[O-])cc1")

                # Check for ethyl ester deprotection
                if product_mol and all(mol for mol in reactant_mols):
                    reactants_have_ethyl_ester = any(
                        mol.HasSubstructMatch(ethyl_ester_pattern) for mol in reactant_mols if mol
                    )
                    product_has_acid = product_mol.HasSubstructMatch(carboxylic_acid_pattern)

                    if reactants_have_ethyl_ester and product_has_acid:
                        print("Detected ethyl ester deprotection")
                        ethyl_ester_deprotection = True

                # Check for p-nitrobenzyl ester protection
                if product_mol and all(mol for mol in reactant_mols):
                    reactants_have_acid = any(
                        mol.HasSubstructMatch(carboxylic_acid_pattern)
                        for mol in reactant_mols
                        if mol
                    )
                    product_has_pnb_ester = product_mol.HasSubstructMatch(pnb_ester_pattern)

                    if reactants_have_acid and product_has_pnb_ester:
                        print("Detected p-nitrobenzyl ester protection")
                        pnb_ester_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both transformations are detected
    return ethyl_ester_deprotection and pnb_ester_protection
