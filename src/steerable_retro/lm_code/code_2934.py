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
    Detects if the synthetic route involves nucleophilic aromatic substitution
    activated by nitro groups.
    """
    found_nitro_activated_snar = False

    def dfs_traverse(node):
        nonlocal found_nitro_activated_snar

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
            aromatic_pattern = Chem.MolFromSmarts("a")

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if (
                        mol
                        and mol.HasSubstructMatch(nitro_pattern)
                        and mol.HasSubstructMatch(aromatic_pattern)
                    ):
                        # Check if this is likely an SNAr reaction
                        # Look for common leaving groups
                        leaving_groups = [
                            Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]"),  # Halides
                            Chem.MolFromSmarts("[#6]-[#16]-[#6]#[#7]"),  # Thiocyanate
                            Chem.MolFromSmarts(
                                "[#6]-[#8]-[#16](=[O])(=[O])[#6]"
                            ),  # Tosylate/mesylate
                        ]

                        for lg_pattern in leaving_groups:
                            if mol.HasSubstructMatch(lg_pattern):
                                prod_mol = Chem.MolFromSmiles(product)
                                if prod_mol and prod_mol.HasSubstructMatch(nitro_pattern):
                                    # If the product still has the nitro group but the leaving group is replaced
                                    if not prod_mol.HasSubstructMatch(lg_pattern) or len(
                                        prod_mol.GetSubstructMatches(lg_pattern)
                                    ) < len(mol.GetSubstructMatches(lg_pattern)):
                                        print("Found nitro-activated SNAr")
                                        found_nitro_activated_snar = True
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_nitro_activated_snar
