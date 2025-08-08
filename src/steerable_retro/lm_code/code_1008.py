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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects if the synthetic route involves methyl ester hydrolysis.
    """
    ester_hydrolyzed = False

    def dfs_traverse(node):
        nonlocal ester_hydrolyzed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction: {rsmi}")

                # Check if this is an ester hydrolysis or related reaction
                is_hydrolysis = (
                    checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                )

                is_esterification = checker.check_reaction(
                    "Esterification of Carboxylic Acids", rsmi
                )

                if is_hydrolysis or is_esterification:
                    reactants_smiles = rsmi.split(">")[0]
                    product_smiles = rsmi.split(">")[-1]

                    # For esterification, we need to swap reactants and products since we're looking at it retrosynthetically
                    if is_esterification:
                        temp = reactants_smiles
                        reactants_smiles = product_smiles
                        product_smiles = temp
                        print(
                            "Swapped reactants and products for esterification (retrosynthetic view)"
                        )

                    # Check if reactants contain methyl ester and product contains carboxylic acid
                    reactants = reactants_smiles.split(".")

                    for r in reactants:
                        if checker.check_fg("Ester", r):
                            print(f"Found ester in reactant: {r}")

                            # Check if product contains carboxylic acid
                            if checker.check_fg("Carboxylic acid", product_smiles):
                                print(f"Found carboxylic acid in product: {product_smiles}")

                                # Check if it's specifically a methyl ester
                                reactant_mol = Chem.MolFromSmiles(r)
                                if reactant_mol:
                                    # Look for methyl group connected to ester oxygen
                                    for atom in reactant_mol.GetAtoms():
                                        if atom.GetSymbol() == "C" and atom.GetDegree() == 1:
                                            for neighbor in atom.GetNeighbors():
                                                if (
                                                    neighbor.GetSymbol() == "O"
                                                    and neighbor.GetDegree() == 2
                                                ):
                                                    for nn in neighbor.GetNeighbors():
                                                        if (
                                                            nn.GetSymbol() == "C"
                                                            and nn.GetIsAromatic() == False
                                                        ):
                                                            for nnn in nn.GetNeighbors():
                                                                if (
                                                                    nnn.GetSymbol() == "O"
                                                                    and nnn.GetDegree() == 1
                                                                ):
                                                                    print(
                                                                        "Confirmed methyl ester hydrolysis"
                                                                    )
                                                                    ester_hydrolyzed = True
                                                                    break

                # Also check for direct esterification in reverse direction
                if not ester_hydrolyzed:
                    reactants_smiles = rsmi.split(">")[0]
                    product_smiles = rsmi.split(">")[-1]

                    if checker.check_fg("Carboxylic acid", reactants_smiles) and checker.check_fg(
                        "Ester", product_smiles
                    ):
                        # Check if it's a methyl ester in the product
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            for atom in product_mol.GetAtoms():
                                if atom.GetSymbol() == "C" and atom.GetDegree() == 1:
                                    for neighbor in atom.GetNeighbors():
                                        if (
                                            neighbor.GetSymbol() == "O"
                                            and neighbor.GetDegree() == 2
                                        ):
                                            for nn in neighbor.GetNeighbors():
                                                if (
                                                    nn.GetSymbol() == "C"
                                                    and nn.GetIsAromatic() == False
                                                ):
                                                    for nnn in nn.GetNeighbors():
                                                        if (
                                                            nnn.GetSymbol() == "O"
                                                            and nnn.GetDegree() == 1
                                                        ):
                                                            print(
                                                                "Confirmed methyl ester formation (reverse of hydrolysis)"
                                                            )
                                                            ester_hydrolyzed = True
                                                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Ester hydrolysis detected: {ester_hydrolyzed}")
    return ester_hydrolyzed
