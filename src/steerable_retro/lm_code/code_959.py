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
    Detects if the synthetic route involves nitrile hydrolysis on a heterocyclic system.
    """
    found_strategy = False

    def dfs_traverse(node):
        nonlocal found_strategy

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product = Chem.MolFromSmiles(product_smiles)

                if all(reactants) and product:
                    # Check for nitrile in reactants
                    nitrile_pattern = Chem.MolFromSmarts("[#6]#[N]")

                    # Check for carboxylic acid in product
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[O])[O;H,-]")

                    # Check for heterocycles
                    imidazole_pattern = Chem.MolFromSmarts("c1ncnn1")
                    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
                    pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

                    # Check if any reactant has nitrile and heterocycle
                    nitrile_on_heterocycle = False
                    for reactant in reactants:
                        has_nitrile = reactant.HasSubstructMatch(nitrile_pattern)
                        has_heterocycle = (
                            reactant.HasSubstructMatch(imidazole_pattern)
                            or reactant.HasSubstructMatch(furan_pattern)
                            or reactant.HasSubstructMatch(pyridine_pattern)
                        )

                        if has_nitrile and has_heterocycle:
                            nitrile_on_heterocycle = True

                    # Check if product has carboxylic acid and heterocycle
                    acid_on_heterocycle = False
                    if product.HasSubstructMatch(carboxylic_acid_pattern) and (
                        product.HasSubstructMatch(imidazole_pattern)
                        or product.HasSubstructMatch(furan_pattern)
                        or product.HasSubstructMatch(pyridine_pattern)
                    ):
                        acid_on_heterocycle = True

                    # If we have nitrile on heterocycle in reactants and acid on heterocycle in product,
                    # this is likely a nitrile hydrolysis on a heterocyclic system
                    if nitrile_on_heterocycle and acid_on_heterocycle:
                        print("Found nitrile hydrolysis on heterocyclic system")
                        found_strategy = True
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_strategy
