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
    This function detects if the synthetic route involves a Suzuki coupling
    for biaryl formation.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have two reactants (typical for Suzuki coupling)
                if len(reactants) == 2:
                    # Try to create molecule objects
                    try:
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                        product_mol = Chem.MolFromSmiles(product)

                        # Check for boronic acid/ester in one reactant
                        boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                        # Check for aryl halide in the other reactant
                        aryl_halide_pattern = Chem.MolFromSmarts("c-[#9,#17,#35,#53]")
                        # Check for biaryl in product
                        biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")

                        has_boronic = any(
                            mol is not None and mol.HasSubstructMatch(boronic_pattern)
                            for mol in reactant_mols
                        )
                        has_aryl_halide = any(
                            mol is not None and mol.HasSubstructMatch(aryl_halide_pattern)
                            for mol in reactant_mols
                        )
                        has_biaryl = product_mol is not None and product_mol.HasSubstructMatch(
                            biaryl_pattern
                        )

                        if has_boronic and has_aryl_halide and has_biaryl:
                            print("Suzuki coupling detected for biaryl formation")
                            suzuki_detected = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return suzuki_detected
