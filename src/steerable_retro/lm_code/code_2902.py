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
    Detects a strategy where an alkyne is first converted to a vinyl sulfide,
    then to a vinyl ether, maintaining the carbon skeleton while changing heteroatom attachments.
    """
    # Initialize tracking variables
    has_alkyne = False
    has_vinyl_sulfide = False
    has_vinyl_ether = False
    reaction_sequence = []

    def dfs_traverse(node):
        nonlocal has_alkyne, has_vinyl_sulfide, has_vinyl_ether, reaction_sequence

        if node["type"] == "mol":
            # Check for functional groups in molecules
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for alkyne
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[C]")):
                    has_alkyne = True
                    print(f"Found alkyne in molecule: {node['smiles']}")

                # Check for vinyl sulfide
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C][S]")):
                    has_vinyl_sulfide = True
                    print(f"Found vinyl sulfide in molecule: {node['smiles']}")

                # Check for vinyl ether
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C][O]")):
                    has_vinyl_ether = True
                    print(f"Found vinyl ether in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            # Analyze reaction transformations
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for thiol addition to alkyne
                if (
                    Chem.MolFromSmiles(reactants).HasSubstructMatch(Chem.MolFromSmarts("[C]#[C]"))
                    and Chem.MolFromSmiles(reactants).HasSubstructMatch(Chem.MolFromSmarts("[SH]"))
                    and Chem.MolFromSmiles(product).HasSubstructMatch(
                        Chem.MolFromSmarts("[C]=[C][S]")
                    )
                ):
                    reaction_sequence.append("alkyne_to_vinyl_sulfide")
                    print(f"Detected alkyne to vinyl sulfide transformation in reaction: {rsmi}")

                # Check for thioether to ether conversion
                if Chem.MolFromSmiles(reactants).HasSubstructMatch(
                    Chem.MolFromSmarts("[C]=[C][S]")
                ) and Chem.MolFromSmiles(product).HasSubstructMatch(
                    Chem.MolFromSmarts("[C]=[C][O]")
                ):
                    reaction_sequence.append("vinyl_sulfide_to_vinyl_ether")
                    print(
                        f"Detected vinyl sulfide to vinyl ether transformation in reaction: {rsmi}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_alkyne
        and has_vinyl_sulfide
        and has_vinyl_ether
        and "alkyne_to_vinyl_sulfide" in reaction_sequence
        and "vinyl_sulfide_to_vinyl_ether" in reaction_sequence
    )

    if strategy_present:
        print("Detected alkyne to vinyl heteroatom strategy")
    else:
        print("Alkyne to vinyl heteroatom strategy not detected")

    return strategy_present
