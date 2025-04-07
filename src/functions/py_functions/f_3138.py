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
    This function detects if the synthetic route involves late-stage introduction
    of an amine functionality via nucleophilic substitution.
    """
    late_stage_amine = False

    def dfs_traverse(node):
        nonlocal late_stage_amine

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 0 or 1)
            if "depth" in node["metadata"] and node["metadata"]["depth"] <= 1:
                try:
                    # Check if any reactant has amine
                    has_amine = False
                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            amine_pattern = Chem.MolFromSmarts("[#7]")
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                                break

                    # Check if any reactant has haloalkyl
                    has_haloalkyl = False
                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            haloalkyl_pattern = Chem.MolFromSmarts("[CH2][Cl,Br,I]")
                            if r_mol.HasSubstructMatch(haloalkyl_pattern):
                                has_haloalkyl = True
                                break

                    # Check if product has C-N bond that wasn't in reactants
                    if has_amine and has_haloalkyl:
                        late_stage_amine = True
                        print(
                            "Detected late-stage amine introduction via nucleophilic substitution at depth",
                            node["metadata"]["depth"],
                        )
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return late_stage_amine
