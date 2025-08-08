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
    Detects a linear synthesis strategy with late-stage Suzuki coupling.
    Key features:
    1. Contains a Suzuki coupling (C-C bond formation with boronic acid/ester)
    2. The Suzuki coupling occurs in the second half of the synthesis
    3. Borylation step precedes the coupling
    """
    suzuki_coupling_found = False
    borylation_found = False
    suzuki_depth = -1
    borylation_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found, borylation_found, suzuki_depth, borylation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(r for r in reactants):
                # Check for borylation (C-Br to C-B transformation)
                has_bromo = any(
                    r.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br]")) for r in reactants
                )
                product_has_boron = product.HasSubstructMatch(Chem.MolFromSmarts("[c]-[B]"))

                if has_bromo and product_has_boron:
                    borylation_found = True
                    borylation_depth = depth
                    print(f"Borylation detected at depth {depth}")

                # Check for Suzuki coupling (biaryl formation with boronic acid/ester reactant)
                has_boron = any(
                    r.HasSubstructMatch(Chem.MolFromSmarts("[c]-[B]")) for r in reactants
                )
                has_halide = any(
                    r.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br,I,Cl]")) for r in reactants
                )

                if has_boron and has_halide:
                    # Further check if a new C-C bond is formed between two aromatic rings
                    suzuki_coupling_found = True
                    suzuki_depth = depth
                    print(f"Suzuki coupling detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Determine if this is a late-stage Suzuki coupling
    # Late stage means in the second half of the synthesis (lower depth values)
    is_late_stage = suzuki_depth <= max_depth / 2

    # Check if borylation precedes Suzuki coupling
    correct_order = borylation_depth > suzuki_depth

    result = suzuki_coupling_found and borylation_found and is_late_stage and correct_order
    print(f"Linear synthesis with late-stage Suzuki: {result}")
    print(
        f"Suzuki depth: {suzuki_depth}, Borylation depth: {borylation_depth}, Max depth: {max_depth}"
    )

    return result
