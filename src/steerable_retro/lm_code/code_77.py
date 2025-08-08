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
    Detects a synthetic strategy where a nitro group is reduced to an amine,
    which is then converted to an amide in the final step, all on an aromatic heterocycle.
    """
    # Initialize tracking variables
    has_nitro_reduction = False
    has_amide_formation = False
    amide_formation_depth = None
    nitro_reduction_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction, has_amide_formation, amide_formation_depth, nitro_reduction_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for nitro reduction
                nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                amine_pattern = Chem.MolFromSmarts("[#7H2]")

                reactant_has_nitro = any(r.HasSubstructMatch(nitro_pattern) for r in reactants if r)
                product_has_amine = product.HasSubstructMatch(amine_pattern) if product else False

                if reactant_has_nitro and product_has_amine:
                    has_nitro_reduction = True
                    nitro_reduction_depth = depth
                    print(f"Found nitro reduction at depth {depth}")

                # Check for amide formation
                amine_pattern = Chem.MolFromSmarts("[#7H2]")
                amide_pattern = Chem.MolFromSmarts("[#7H]-[#6](=[#8])")

                reactant_has_amine = any(r.HasSubstructMatch(amine_pattern) for r in reactants if r)
                product_has_amide = product.HasSubstructMatch(amide_pattern) if product else False

                if reactant_has_amine and product_has_amide:
                    has_amide_formation = True
                    amide_formation_depth = depth
                    print(f"Found amide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # The strategy requires both nitro reduction and amide formation,
    # with amide formation occurring at a lower depth (later in synthesis)
    strategy_present = (
        has_nitro_reduction
        and has_amide_formation
        and (amide_formation_depth is not None)
        and (nitro_reduction_depth is not None)
        and (amide_formation_depth < nitro_reduction_depth)
    )

    print(f"Nitro-to-amide strategy detected: {strategy_present}")
    return strategy_present
