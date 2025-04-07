#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the synthetic route uses chlorine as a leaving group in nucleophilic substitutions.
    """
    chlorine_leaving_group_used = False

    def dfs_traverse(node):
        nonlocal chlorine_leaving_group_used

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for C-Cl bond in reactants
                c_cl_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
                reactant_has_c_cl = False
                for r_smiles in reactants_smiles:
                    try:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(c_cl_pattern):
                            reactant_has_c_cl = True
                            break
                    except:
                        continue

                # Check if C-Cl bond is broken and C-N bond is formed (nucleophilic substitution)
                if reactant_has_c_cl:
                    try:
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        # Check if product has fewer C-Cl bonds than reactants
                        # and if it has a new C-N bond
                        c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                        if product_mol:
                            # Simple heuristic: if product has C-N bond and reactant had C-Cl
                            # that's now gone, it's likely a substitution
                            if product_mol.HasSubstructMatch(c_n_pattern):
                                print(
                                    "Chlorine leaving group in nucleophilic substitution detected"
                                )
                                chlorine_leaving_group_used = True
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return chlorine_leaving_group_used
