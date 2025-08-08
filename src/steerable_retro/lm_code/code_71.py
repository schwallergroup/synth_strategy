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
    This function detects if the synthetic route involves a late-stage (depth 0 or 1)
    aromatic C-N bond formation, typically via nucleophilic aromatic substitution.
    """
    late_stage_cn_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cn_formation

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if an amine is present in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;!$(N~[!#6]);!$(N~[#6]~[#7,#8,#16])]")
            for reactant in reactants_smiles.split("."):
                try:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(amine_pattern):
                        # Check if product has new C-N bond where C is aromatic
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            aromatic_cn_pattern = Chem.MolFromSmarts("c-[#7]")
                            if product_mol.HasSubstructMatch(aromatic_cn_pattern):
                                print(
                                    f"Found late-stage aromatic C-N bond formation at depth {depth}"
                                )
                                late_stage_cn_formation = True
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_cn_formation
