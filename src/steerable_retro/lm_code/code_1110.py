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
    This function detects if the synthetic route involves O-methylation of phenol.
    """
    phenol_pattern = Chem.MolFromSmarts("c[OH]")
    methoxy_pattern = Chem.MolFromSmarts("cO[CH3]")
    o_methylation_detected = False

    def dfs_traverse(node):
        nonlocal o_methylation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant contains phenol
            reactant_has_phenol = False
            for r_smiles in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and r_mol.HasSubstructMatch(phenol_pattern):
                        reactant_has_phenol = True
                        break
                except:
                    continue

            # Check if product contains methoxy group
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                product_has_methoxy = p_mol and p_mol.HasSubstructMatch(methoxy_pattern)
            except:
                product_has_methoxy = False

            # If reactant has phenol and product has methoxy, it's an O-methylation
            if reactant_has_phenol and product_has_methoxy:
                print(f"O-methylation detected in reaction: {rsmi}")
                o_methylation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return o_methylation_detected
