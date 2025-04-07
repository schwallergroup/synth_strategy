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
    Detects a synthetic strategy involving multiple SNAr reactions,
    particularly focusing on late-stage SNAr on pyrimidine core.
    """
    snar_reactions_count = 0
    late_stage_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_reactions_count, late_stage_snar

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr pattern: halogen on aromatic being replaced by O or N
                # Look for patterns where Cl/F/Br on aromatic is replaced by O/N
                product_mol = Chem.MolFromSmiles(product)

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and product_mol:
                        # Check for aromatic halogen in reactant
                        ar_hal_pattern = Chem.MolFromSmarts("c-[Cl,F,Br]")
                        if reactant_mol.HasSubstructMatch(ar_hal_pattern):
                            # Check for O/N connection in product where halogen was
                            ar_o_pattern = Chem.MolFromSmarts("c-[O,N]")
                            if product_mol.HasSubstructMatch(ar_o_pattern):
                                # Check if it's an electron-deficient ring (pyrimidine or with nitro)
                                el_def_pattern1 = Chem.MolFromSmarts(
                                    "c1ncncc1"
                                )  # pyrimidine
                                el_def_pattern2 = Chem.MolFromSmarts(
                                    "c-[N+](=[O])-[O-]"
                                )  # nitro

                                if reactant_mol.HasSubstructMatch(
                                    el_def_pattern1
                                ) or reactant_mol.HasSubstructMatch(el_def_pattern2):
                                    snar_reactions_count += 1

                                    # Check if it's late stage (depth 0 or 1)
                                    if depth <= 1:
                                        late_stage_snar = True
                                        print(
                                            f"Late-stage SNAr detected at depth {depth}"
                                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy criteria: at least 2 SNAr reactions with at least one in late stage
    result = snar_reactions_count >= 2 and late_stage_snar
    print(
        f"Sequential SNAr strategy detected: {result} (count: {snar_reactions_count}, late stage: {late_stage_snar})"
    )
    return result
