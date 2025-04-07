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
    Detects a synthetic route that includes amide coupling between a carboxylic acid and an amine,
    where at least one component is a heterocycle.
    """
    has_heterocycle_amide_coupling = False

    def dfs_traverse(node):
        nonlocal has_heterocycle_amide_coupling

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            amide_pattern = Chem.MolFromSmarts("NC(=O)")

            # Heterocycle patterns
            pyrrolo_pyridine_pattern = Chem.MolFromSmarts("c1cnc2[nH]ccc2c1")
            imidazo_pyrimidine_pattern = Chem.MolFromSmarts("c1ncn2ccnc2c1")

            has_acid = any(
                r and r.HasSubstructMatch(carboxylic_acid_pattern) for r in reactant_mols
            )
            has_amine = any(r and r.HasSubstructMatch(amine_pattern) for r in reactant_mols)
            has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

            has_heterocycle = any(
                r
                and (
                    r.HasSubstructMatch(pyrrolo_pyridine_pattern)
                    or r.HasSubstructMatch(imidazo_pyrimidine_pattern)
                )
                for r in reactant_mols
            )

            if has_acid and has_amine and has_amide and has_heterocycle:
                has_heterocycle_amide_coupling = True
                print("Found heterocycle amide coupling")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Heterocycle amide coupling strategy: {has_heterocycle_amide_coupling}")
    return has_heterocycle_amide_coupling
