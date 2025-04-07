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
    This function detects if the synthetic route involves a nitro reduction followed by amide formation.
    """
    # Track molecules through the synthesis
    molecule_sequence = []
    nitro_reduction_step = -1
    amide_formation_step = -1

    def dfs_traverse(node, step=0):
        nonlocal nitro_reduction_step, amide_formation_step

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction
                nitro_pattern = Chem.MolFromSmarts("[#6]c[N+](=[O])[O-]")
                amine_pattern = Chem.MolFromSmarts("[#6]c[NH2]")

                nitro_in_reactants = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(nitro_pattern):
                            nitro_in_reactants = True
                            break
                    except:
                        continue

                if nitro_in_reactants:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(amine_pattern):
                            print("Nitro reduction detected at step", step)
                            nitro_reduction_step = step
                    except:
                        pass

                # Check for amide formation
                acyl_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")
                amide_pattern = Chem.MolFromSmarts("C(=O)N")

                acyl_chloride_in_reactants = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(acyl_chloride_pattern):
                            acyl_chloride_in_reactants = True
                            break
                    except:
                        continue

                if acyl_chloride_in_reactants:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(amide_pattern):
                            print("Amide formation detected at step", step)
                            amide_formation_step = step
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child, step + 1)

    dfs_traverse(route)

    # Check if nitro reduction occurs before amide formation
    if nitro_reduction_step != -1 and amide_formation_step != -1:
        if (
            nitro_reduction_step > amide_formation_step
        ):  # Remember: higher step number = earlier in synthesis
            print("Sequence detected: nitro reduction followed by amide formation")
            return True

    return False
