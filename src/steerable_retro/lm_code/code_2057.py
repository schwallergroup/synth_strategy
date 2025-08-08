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
    This function detects a sequence of functional group interconversions:
    ester → amide → amine → amide
    """
    # Define SMARTS patterns for functional groups
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    amide_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[NX3]")
    amine_pattern = Chem.MolFromSmarts("[#6][NX3]")

    # Track the sequence of functional group transformations
    transformations = []

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Determine functional groups in reactants and products
                product_has_ester = product_mol and product_mol.HasSubstructMatch(ester_pattern)
                product_has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)
                product_has_amine = product_mol and product_mol.HasSubstructMatch(amine_pattern)

                reactants_have_ester = any(
                    r and r.HasSubstructMatch(ester_pattern) for r in reactant_mols if r
                )
                reactants_have_amide = any(
                    r and r.HasSubstructMatch(amide_pattern) for r in reactant_mols if r
                )
                reactants_have_amine = any(
                    r and r.HasSubstructMatch(amine_pattern) for r in reactant_mols if r
                )

                # Determine transformation type
                if reactants_have_ester and product_has_amide and not product_has_ester:
                    transformations.append("ester_to_amide")
                    print(
                        f"Found ester to amide transformation at {node['metadata'].get('ID', '')}"
                    )
                elif reactants_have_amide and product_has_amine and not product_has_amide:
                    transformations.append("amide_to_amine")
                    print(
                        f"Found amide to amine transformation at {node['metadata'].get('ID', '')}"
                    )
                elif reactants_have_amine and product_has_amide and not reactants_have_amide:
                    transformations.append("amine_to_amide")
                    print(
                        f"Found amine to amide transformation at {node['metadata'].get('ID', '')}"
                    )
            except:
                print("Error processing SMILES in functional group detection")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the desired sequence
    sequence_str = "_".join(transformations)
    return (
        "ester_to_amide" in sequence_str
        and "amide_to_amine" in sequence_str
        and "amine_to_amide" in sequence_str
    )
