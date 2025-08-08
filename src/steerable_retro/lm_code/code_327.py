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
    This function detects a synthetic strategy involving ester hydrolysis followed by amide coupling.
    It looks for a sequence where an ester is hydrolyzed to a carboxylic acid, which is then used
    in an amide coupling reaction.
    """
    # Track if we found the pattern
    found_ester_hydrolysis = False
    found_amide_coupling = False

    # Track the depth of each reaction
    ester_hydrolysis_depth = -1
    amide_coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_hydrolysis, found_amide_coupling
        nonlocal ester_hydrolysis_depth, amide_coupling_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester hydrolysis
            if not found_ester_hydrolysis:
                # Look for ester pattern in reactants
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
                # Look for carboxylic acid pattern in product
                acid_pattern = Chem.MolFromSmarts("[C](=[O])[O]")

                reactant_has_ester = False
                for r_smiles in reactants_smiles:
                    try:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(ester_pattern):
                            reactant_has_ester = True
                            break
                    except:
                        continue

                product_has_acid = False
                try:
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    if p_mol and p_mol.HasSubstructMatch(acid_pattern):
                        product_has_acid = True
                except:
                    pass

                if reactant_has_ester and product_has_acid:
                    found_ester_hydrolysis = True
                    ester_hydrolysis_depth = depth
                    print(f"Found ester hydrolysis at depth {depth}")

            # Check for amide coupling
            if not found_amide_coupling:
                # Look for carboxylic acid pattern in reactants
                acid_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
                # Look for amine pattern in reactants
                amine_pattern = Chem.MolFromSmarts("[N]")
                # Look for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[N][C](=[O])")

                reactant_has_acid = False
                reactant_has_amine = False

                for r_smiles in reactants_smiles:
                    try:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol:
                            if r_mol.HasSubstructMatch(acid_pattern):
                                reactant_has_acid = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                reactant_has_amine = True
                    except:
                        continue

                product_has_amide = False
                try:
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    if p_mol and p_mol.HasSubstructMatch(amide_pattern):
                        product_has_amide = True
                except:
                    pass

                if reactant_has_acid and reactant_has_amine and product_has_amide:
                    found_amide_coupling = True
                    amide_coupling_depth = depth
                    print(f"Found amide coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the pattern and in the correct order
    if found_ester_hydrolysis and found_amide_coupling:
        # Ester hydrolysis should happen before amide coupling (higher depth)
        if ester_hydrolysis_depth > amide_coupling_depth:
            print("Found ester hydrolysis followed by amide coupling strategy")
            return True

    return False
