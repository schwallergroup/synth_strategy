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
    Detects if the synthetic route involves N-arylation of an indazole core.
    """
    n_arylation_found = False

    def dfs_traverse(node):
        nonlocal n_arylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Pattern for N-arylated indazole
                n_arylated_indazole = Chem.MolFromSmarts(
                    "[c]1[c][c][c]([c][c]1)-[n]1[n][c][c]2[c]1[c][c][c][c]2"
                )

                # Pattern for NH-indazole
                nh_indazole = Chem.MolFromSmarts("[nH]1[n][c][c]2[c]1[c][c][c][c]2")

                # Pattern for aryl boronic acid or halide
                aryl_coupling_partner = Chem.MolFromSmarts("[c][B,Br,I,Cl,Sn]")

                if (
                    product_mol
                    and n_arylated_indazole
                    and product_mol.HasSubstructMatch(n_arylated_indazole)
                ):
                    # Check if reactants contain NH-indazole and aryl coupling partner
                    has_nh_indazole = False
                    has_aryl_partner = False

                    for r_mol in reactant_mols:
                        if r_mol:
                            if nh_indazole and r_mol.HasSubstructMatch(nh_indazole):
                                has_nh_indazole = True
                            if aryl_coupling_partner and r_mol.HasSubstructMatch(
                                aryl_coupling_partner
                            ):
                                has_aryl_partner = True

                    if has_nh_indazole and has_aryl_partner:
                        n_arylation_found = True
                        print(f"Found N-arylation of indazole: {rsmi}")
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return n_arylation_found
