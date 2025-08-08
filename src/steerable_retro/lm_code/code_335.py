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
    Detects a synthesis strategy where two different nitrogen heterocycles
    are connected using a bifunctional alkyl chain linker through sequential N-alkylation reactions.
    """
    # Track if we found the required patterns
    found_phthalimide_alkylation = False
    found_morpholine_alkylation = False
    found_dibromoalkane = False

    # SMARTS patterns
    phthalimide_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)NC(=O)")
    morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")
    dibromoalkane_pattern = Chem.MolFromSmarts("[Br][#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6][Br]")

    def dfs_traverse(node):
        nonlocal found_phthalimide_alkylation, found_morpholine_alkylation, found_dibromoalkane

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    # Check for dibromoalkane in reactants
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(dibromoalkane_pattern):
                            found_dibromoalkane = True
                            print("Found dibromoalkane reactant")

                    # Check for phthalimide alkylation
                    # Look for a reaction where phthalimide is in reactants and
                    # alkylated phthalimide is in product
                    phthalimide_in_reactants = any(
                        r_mol and r_mol.HasSubstructMatch(phthalimide_pattern)
                        for r_mol in reactant_mols
                    )
                    phthalimide_in_product = product_mol and product_mol.HasSubstructMatch(
                        phthalimide_pattern
                    )

                    if phthalimide_in_reactants and phthalimide_in_product:
                        # Check if this is an alkylation (N-C bond formation)
                        # This is a simplification - in a real implementation, you would need to
                        # analyze the reaction more carefully
                        if any("Br" in r for r in reactants):
                            found_phthalimide_alkylation = True
                            print("Found phthalimide alkylation")

                    # Check for morpholine alkylation
                    morpholine_in_reactants = any(
                        r_mol and r_mol.HasSubstructMatch(morpholine_pattern)
                        for r_mol in reactant_mols
                    )
                    morpholine_in_product = product_mol and product_mol.HasSubstructMatch(
                        morpholine_pattern
                    )

                    if morpholine_in_reactants and morpholine_in_product:
                        # Check if this is an alkylation (N-C bond formation)
                        if any("Br" in r for r in reactants):
                            found_morpholine_alkylation = True
                            print("Found morpholine alkylation")

                except:
                    print("Error processing SMILES in reaction")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if all required patterns were found
    strategy_detected = (
        found_phthalimide_alkylation and found_morpholine_alkylation and found_dibromoalkane
    )
    print(f"Sequential heterocycle connection strategy detected: {strategy_detected}")
    return strategy_detected
