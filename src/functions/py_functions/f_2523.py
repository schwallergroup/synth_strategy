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
    Detects if the synthesis uses a late-stage fragment coupling strategy
    with N-acylation as the key bond-forming reaction.
    """
    found_late_acylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_acylation

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an acylation reaction
                acyl_pattern = Chem.MolFromSmarts("[#7:1][C:2](=[O:3])[#6:4]")
                prod_mol = Chem.MolFromSmiles(product)

                if prod_mol and prod_mol.HasSubstructMatch(acyl_pattern):
                    # Check if reactants include an acid chloride or similar acylating agent
                    acyl_agent_pattern = Chem.MolFromSmarts("[Cl,Br,I,OH][C](=[O])[#6]")
                    amine_pattern = Chem.MolFromSmarts("[#7;!$(N~[!#6])]")

                    has_acyl_agent = False
                    has_amine = False

                    for reactant in reactants:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            if r_mol.HasSubstructMatch(acyl_agent_pattern):
                                has_acyl_agent = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    if has_acyl_agent and has_amine:
                        # Check if both fragments are complex (have multiple atoms)
                        complex_fragments = 0
                        for reactant in reactants:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if (
                                r_mol and r_mol.GetNumAtoms() > 5
                            ):  # Arbitrary threshold for "complex"
                                complex_fragments += 1

                        if complex_fragments >= 2:
                            found_late_acylation = True
                            print(
                                f"Found late-stage acylation coupling at depth {depth}"
                            )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_late_acylation
