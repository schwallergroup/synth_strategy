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
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects sulfilimine (S=N) formation from thioether.
    Note: Despite the function name, this actually detects the oxidation of thioether to sulfilimine.
    """
    sulfilimine_formation_detected = False

    def dfs_traverse(node):
        nonlocal sulfilimine_formation_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product_part = rsmi.split(">")[-1]
            products = product_part.split(".")

            print(f"Analyzing reaction: {rsmi}")

            # Check for thioether pattern in reactants
            thioether_reactant = None
            for reactant in reactants:
                if checker.check_fg("Monosulfide", reactant):
                    thioether_reactant = reactant
                    print(f"Found thioether in reactant: {reactant}")
                    break

            # If no thioether found, continue to next node
            if not thioether_reactant:
                print("No thioether found in reactants")
                for child in node.get("children", []):
                    dfs_traverse(child)
                return

            # Check for sulfilimine pattern in products
            sulfilimine_product = None
            for prod in products:
                # Look for S=N pattern (sulfilimine)
                mol = Chem.MolFromSmiles(prod)
                if mol:
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "S":
                            for bond in atom.GetBonds():
                                other_atom = bond.GetOtherAtom(atom)
                                if (
                                    other_atom.GetSymbol() == "N"
                                    and bond.GetBondType() == Chem.BondType.DOUBLE
                                ):
                                    sulfilimine_product = prod
                                    print(f"Found sulfilimine in product: {prod}")
                                    break
                        if sulfilimine_product:
                            break

            # If no sulfilimine found, continue to next node
            if not sulfilimine_product:
                print("No sulfilimine found in products")
                for child in node.get("children", []):
                    dfs_traverse(child)
                return

            # Verify the S atom in thioether maps to S in sulfilimine using atom mapping
            if thioether_reactant and sulfilimine_product:
                thioether_mol = Chem.MolFromSmiles(thioether_reactant)
                sulfilimine_mol = Chem.MolFromSmiles(sulfilimine_product)

                if thioether_mol and sulfilimine_mol:
                    # Look for mapped sulfur atoms
                    thioether_s_maps = []
                    for atom in thioether_mol.GetAtoms():
                        if atom.GetSymbol() == "S" and atom.GetAtomMapNum() > 0:
                            thioether_s_maps.append(atom.GetAtomMapNum())

                    sulfilimine_s_maps = []
                    for atom in sulfilimine_mol.GetAtoms():
                        if atom.GetSymbol() == "S" and atom.GetAtomMapNum() > 0:
                            sulfilimine_s_maps.append(atom.GetAtomMapNum())

                    # Check if any sulfur atom is mapped between reactant and product
                    common_maps = set(thioether_s_maps).intersection(
                        set(sulfilimine_s_maps)
                    )
                    if common_maps:
                        print(f"Found common mapped S atoms: {common_maps}")
                        sulfilimine_formation_detected = True
                        print("Sulfilimine formation confirmed!")
                    else:
                        # If no atom mapping or common maps found, assume it's a valid transformation
                        # This is a fallback for cases where atom mapping is not available
                        if not thioether_s_maps or not sulfilimine_s_maps:
                            print(
                                "No atom mapping found, assuming valid transformation based on functional groups"
                            )
                            sulfilimine_formation_detected = True
                            print("Sulfilimine formation assumed!")

        for child in node.get("children", []):
            if not sulfilimine_formation_detected:  # Only continue if not yet found
                dfs_traverse(child)

    dfs_traverse(route)
    return sulfilimine_formation_detected
