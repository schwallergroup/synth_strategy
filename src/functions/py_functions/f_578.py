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
    This function detects a synthetic strategy involving fused heterocycle construction
    followed by late-stage SNAr functionalization.
    """
    # Track if we found the key features
    has_fused_heterocycle_formation = False
    has_chlorination = False
    has_snar = False
    has_urea_formation = False

    # Track the depth of reactions
    max_depth = 0
    chlorination_depth = None
    snar_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_fused_heterocycle_formation, has_chlorination, has_snar
        nonlocal has_urea_formation, max_depth, chlorination_depth, snar_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for fused heterocycle formation (oxazole + pyrimidine)
                oxazolo_pyrimidine_pattern = Chem.MolFromSmarts("c1nc2occc2cn1")
                if (
                    product
                    and oxazolo_pyrimidine_pattern
                    and product.HasSubstructMatch(oxazolo_pyrimidine_pattern)
                ):
                    # Check if this is a ring formation by comparing ring counts
                    product_rings = product.GetRingInfo().NumRings()
                    reactant_rings_total = sum(
                        r.GetRingInfo().NumRings() for r in reactants if r
                    )

                    if product_rings > reactant_rings_total:
                        has_fused_heterocycle_formation = True
                        print(f"Found fused heterocycle formation at depth {depth}")

                # Check for chlorination (carbonyl to chloride)
                carbonyl_pattern = Chem.MolFromSmarts("[#8]=[#6]")
                chloride_pattern = Chem.MolFromSmarts("[Cl]-[#6]")

                reactants_have_carbonyl = any(
                    r and r.HasSubstructMatch(carbonyl_pattern) for r in reactants if r
                )
                product_has_chloride = product and product.HasSubstructMatch(
                    chloride_pattern
                )

                if reactants_have_carbonyl and product_has_chloride:
                    has_chlorination = True
                    chlorination_depth = depth
                    print(f"Found chlorination at depth {depth}")

                # Check for SNAr (chloride displacement by amine)
                amine_pattern = Chem.MolFromSmarts("[#7;H2]")
                c_n_pattern = Chem.MolFromSmarts("[#6]-[#7;!$(N=*)]")

                reactants_have_chloride = any(
                    r and r.HasSubstructMatch(chloride_pattern) for r in reactants if r
                )
                reactants_have_amine = any(
                    r and r.HasSubstructMatch(amine_pattern) for r in reactants if r
                )
                product_has_c_n = product and product.HasSubstructMatch(c_n_pattern)

                if reactants_have_chloride and reactants_have_amine and product_has_c_n:
                    has_snar = True
                    snar_depth = depth
                    print(f"Found SNAr at depth {depth}")

                # Check for urea formation
                urea_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#7]")
                product_has_urea = product and product.HasSubstructMatch(urea_pattern)

                if product_has_urea and not any(
                    r and r.HasSubstructMatch(urea_pattern) for r in reactants if r
                ):
                    has_urea_formation = True
                    print(f"Found urea formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # The strategy requires fused heterocycle formation, chlorination, and SNAr
    # The SNAr should be at a lower depth (later in synthesis) than chlorination
    strategy_present = (
        has_fused_heterocycle_formation
        and has_chlorination
        and has_snar
        and (
            snar_depth is not None
            and chlorination_depth is not None
            and snar_depth < chlorination_depth
        )
    )

    print(
        f"Fused heterocycle with late-stage SNAr strategy detected: {strategy_present}"
    )
    return strategy_present
