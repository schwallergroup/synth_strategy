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
    This function detects a synthesis strategy where a complex scaffold is preserved
    while a single position undergoes multiple functional group interconversions.
    """
    # Track functional group changes
    fg_sequence = []
    scaffold_preserved = True

    def dfs_traverse(node):
        nonlocal fg_sequence, scaffold_preserved

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for functional group transformations
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Detect functional groups
                    alcohol_patt = Chem.MolFromSmarts("[#6]-[#8;H1]")
                    aldehyde_patt = Chem.MolFromSmarts("[#6;H1]=O")
                    chloride_patt = Chem.MolFromSmarts("[#6]-[Cl]")
                    nitrile_patt = Chem.MolFromSmarts("[#6]-[#6]#[#7]")
                    carboxylic_patt = Chem.MolFromSmarts("[#6]-[#6](=O)-[#8;H1]")
                    amide_patt = Chem.MolFromSmarts("[#6]-[#6](=O)-[#7]")

                    # Check which functional groups are present
                    r_has_alcohol = reactant_mol.HasSubstructMatch(alcohol_patt)
                    r_has_aldehyde = reactant_mol.HasSubstructMatch(aldehyde_patt)
                    r_has_chloride = reactant_mol.HasSubstructMatch(chloride_patt)
                    r_has_nitrile = reactant_mol.HasSubstructMatch(nitrile_patt)
                    r_has_carboxylic = reactant_mol.HasSubstructMatch(carboxylic_patt)
                    r_has_amide = reactant_mol.HasSubstructMatch(amide_patt)

                    p_has_alcohol = product_mol.HasSubstructMatch(alcohol_patt)
                    p_has_aldehyde = product_mol.HasSubstructMatch(aldehyde_patt)
                    p_has_chloride = product_mol.HasSubstructMatch(chloride_patt)
                    p_has_nitrile = product_mol.HasSubstructMatch(nitrile_patt)
                    p_has_carboxylic = product_mol.HasSubstructMatch(carboxylic_patt)
                    p_has_amide = product_mol.HasSubstructMatch(amide_patt)

                    # Determine the transformation
                    if r_has_alcohol and p_has_aldehyde:
                        fg_sequence.append("alcohol_to_aldehyde")
                    elif r_has_aldehyde and p_has_alcohol:
                        fg_sequence.append("aldehyde_to_alcohol")
                    elif r_has_alcohol and p_has_chloride:
                        fg_sequence.append("alcohol_to_chloride")
                    elif r_has_chloride and p_has_nitrile:
                        fg_sequence.append("chloride_to_nitrile")
                    elif r_has_nitrile and p_has_carboxylic:
                        fg_sequence.append("nitrile_to_carboxylic")
                    elif r_has_carboxylic and p_has_amide:
                        fg_sequence.append("carboxylic_to_amide")

                    # Check if core scaffold is preserved
                    # Define core scaffold SMARTS (indazole + cyclohexyl + methylsulfonylphenyl)
                    scaffold_smarts = "[#6]1:[#6]:[#7]:[#6]2:[#7]:[#6]:[#6]:[#6]:2:[#6]:1"
                    scaffold_patt = Chem.MolFromSmarts(scaffold_smarts)

                    if not (
                        reactant_mol.HasSubstructMatch(scaffold_patt)
                        and product_mol.HasSubstructMatch(scaffold_patt)
                    ):
                        scaffold_preserved = False

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if this matches our strategy
    has_alcohol_aldehyde_bidirectional = (
        "alcohol_to_aldehyde" in fg_sequence and "aldehyde_to_alcohol" in fg_sequence
    )

    has_nitrile_acid_amide_sequence = all(
        x in fg_sequence for x in ["nitrile_to_carboxylic", "carboxylic_to_amide"]
    )

    print(f"FG Sequence: {fg_sequence}")
    print(f"Scaffold preserved: {scaffold_preserved}")
    print(f"Has alcohol-aldehyde bidirectional: {has_alcohol_aldehyde_bidirectional}")
    print(f"Has nitrile-acid-amide sequence: {has_nitrile_acid_amide_sequence}")

    return (
        scaffold_preserved
        and len(fg_sequence) >= 5  # Multiple FG transformations
        and has_alcohol_aldehyde_bidirectional
        and has_nitrile_acid_amide_sequence
    )
