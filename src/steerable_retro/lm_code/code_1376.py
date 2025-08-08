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
    Detects if the synthesis route involves elaboration of heterocycles through
    sequential functionalization (e.g., halogenation followed by coupling).
    """
    # Track heterocycle modifications
    heterocycle_mods = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Define heterocycle patterns
                thiazole_pattern = Chem.MolFromSmarts("[#6]1[#16][#6][#7][#6]1")
                pyrimidine_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#6][#6][#7]1")

                # Define modification patterns
                halogenation_pattern = Chem.MolFromSmarts("[#6]-[Br,Cl,I]")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Check for heterocycle halogenation
                    has_heterocycle_reactant = any(
                        r.HasSubstructMatch(thiazole_pattern)
                        or r.HasSubstructMatch(pyrimidine_pattern)
                        for r in reactants
                        if r
                    )

                    has_halogenated_product = product.HasSubstructMatch(halogenation_pattern)

                    if has_heterocycle_reactant and has_halogenated_product:
                        heterocycle_mods.append(("halogenation", depth))
                        print(f"Heterocycle halogenation detected at depth {depth}")

                    # Check for coupling on halogenated heterocycle
                    has_halogenated_reactant = any(
                        r.HasSubstructMatch(halogenation_pattern) for r in reactants if r
                    )

                    if has_halogenated_reactant and has_heterocycle_reactant:
                        # Check if halogen is removed in product (indicating coupling)
                        halogen_count_reactants = sum(
                            len(r.GetSubstructMatches(halogenation_pattern)) for r in reactants if r
                        )
                        halogen_count_product = len(
                            product.GetSubstructMatches(halogenation_pattern)
                        )

                        if halogen_count_product < halogen_count_reactants:
                            heterocycle_mods.append(("coupling", depth))
                            print(f"Heterocycle coupling detected at depth {depth}")
            except:
                print(f"Error processing reaction SMILES at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have a sequence of halogenation followed by coupling
    halogenation_depths = [
        depth for mod_type, depth in heterocycle_mods if mod_type == "halogenation"
    ]
    coupling_depths = [depth for mod_type, depth in heterocycle_mods if mod_type == "coupling"]

    if halogenation_depths and coupling_depths:
        # Check if halogenation occurs before coupling in the synthesis direction
        # (which means higher depth in retrosynthesis)
        for h_depth in halogenation_depths:
            for c_depth in coupling_depths:
                if h_depth > c_depth:  # In retrosynthesis, higher depth = earlier in synthesis
                    print(
                        "Detected heterocycle elaboration strategy: halogenation followed by coupling"
                    )
                    return True

    return False
