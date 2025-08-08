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
    This function detects late-stage introduction of an amine via nucleophilic aromatic substitution.
    """
    late_stage_amine_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amine_found

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for chloro-aromatic in reactants
                chloro_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
                amine_pattern = Chem.MolFromSmarts(
                    "[#7H][#6][#6][#7]"
                )  # Pattern for aliphatic amine

                # Check if one reactant has chloro and another has amine
                has_chloro = False
                has_amine = False

                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            if r_mol.HasSubstructMatch(chloro_pattern):
                                has_chloro = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                    except:
                        continue

                # Check if product has new C-N bond where chlorine was
                if has_chloro and has_amine:
                    try:
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol:
                            # This is a simplification - in reality would need reaction mapping
                            if p_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#7][#6][#6][#7]")):
                                print(
                                    "Detected late-stage amine introduction via nucleophilic aromatic substitution"
                                )
                                late_stage_amine_found = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return late_stage_amine_found
