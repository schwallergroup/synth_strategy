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
    Detects if the synthetic route involves multiple amide formation reactions.
    """
    amide_formation_count = 0

    def dfs_traverse(node):
        nonlocal amide_formation_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            acyl_pattern = Chem.MolFromSmarts("C(=O)[Cl,O]")
            amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")

            try:
                # Check if reactants contain amine and acyl group
                reactant_has_amine = False
                reactant_has_acyl = False

                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        if r_mol.HasSubstructMatch(amine_pattern):
                            reactant_has_amine = True
                        if r_mol.HasSubstructMatch(acyl_pattern):
                            reactant_has_acyl = True

                # Check if product contains amide
                p_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_amine
                    and reactant_has_acyl
                    and p_mol
                    and p_mol.HasSubstructMatch(amide_pattern)
                ):
                    print(f"Amide formation detected: {rsmi}")
                    amide_formation_count += 1
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return amide_formation_count >= 2
