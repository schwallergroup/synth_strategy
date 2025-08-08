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
    This function detects if the synthesis route involves the transformation of a nitrile group
    to a lactone structure, which is a key functional group interconversion.
    """
    nitrile_to_lactone = False

    # Helper function to detect lactone structures
    def has_lactone_structure(mol_smiles):
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Look for carbonyl carbon that's in a ring and connected to ring oxygen
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "C" and atom.IsInRing():
                has_carbonyl = False
                has_ring_ether = False

                for bond in atom.GetBonds():
                    other_atom = bond.GetOtherAtom(atom)
                    # Check for carbonyl
                    if bond.GetBondType() == Chem.BondType.DOUBLE and other_atom.GetSymbol() == "O":
                        has_carbonyl = True
                    # Check for ether oxygen in ring
                    if (
                        bond.GetBondType() == Chem.BondType.SINGLE
                        and other_atom.GetSymbol() == "O"
                        and other_atom.IsInRing()
                    ):
                        has_ring_ether = True

                if has_carbonyl and has_ring_ether:
                    return True
        return False

    def dfs_traverse(node):
        nonlocal nitrile_to_lactone

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for nitrile in reactants
                    nitrile_found = False
                    nitrile_atom_maps = set()

                    for reactant in reactants:
                        if checker.check_fg("Nitrile", reactant):
                            print(f"Nitrile found in reactant: {reactant}")
                            nitrile_found = True

                            # Get atom mapping numbers for the nitrile carbon
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            for atom in reactant_mol.GetAtoms():
                                # Nitrile carbon is connected to nitrogen with triple bond
                                if atom.GetSymbol() == "C" and atom.GetIsAromatic() == False:
                                    for bond in atom.GetBonds():
                                        other_atom = bond.GetOtherAtom(atom)
                                        if (
                                            other_atom.GetSymbol() == "N"
                                            and bond.GetBondType() == Chem.BondType.TRIPLE
                                        ):
                                            map_num = (
                                                atom.GetProp("molAtomMapNumber")
                                                if atom.HasProp("molAtomMapNumber")
                                                else None
                                            )
                                            if map_num:
                                                nitrile_atom_maps.add(map_num)
                                                print(f"Nitrile carbon atom map: {map_num}")

                    # Check for lactone in product
                    if nitrile_found:
                        product_mol = Chem.MolFromSmiles(product)

                        # Direct check for lactone structure
                        is_lactone = has_lactone_structure(product)
                        print(f"Product has lactone structure: {is_lactone}")

                        # Check for various ring structures that could form lactones
                        ring_types = ["oxolane", "oxane", "dioxane", "dioxolane", "oxetane"]
                        has_ring = any(checker.check_ring(ring, product) for ring in ring_types)
                        has_ester = checker.check_fg("Ester", product)

                        print(f"Product has relevant ring: {has_ring}, ester: {has_ester}")

                        if is_lactone:
                            print(f"Lactone found in product: {product}")

                            # Verify the nitrile carbon is part of the ring structure
                            ring_carbon_maps = set()

                            # Find ring carbons in the product
                            for atom in product_mol.GetAtoms():
                                if atom.GetSymbol() == "C" and atom.HasProp("molAtomMapNumber"):
                                    if atom.IsInRing():
                                        map_num = atom.GetProp("molAtomMapNumber")
                                        ring_carbon_maps.add(map_num)

                            print(f"Ring carbon maps: {ring_carbon_maps}")
                            print(f"Nitrile carbon maps: {nitrile_atom_maps}")

                            # Check if any nitrile carbon maps are in the ring
                            if (
                                nitrile_atom_maps
                                and ring_carbon_maps
                                and nitrile_atom_maps.intersection(ring_carbon_maps)
                            ):
                                print(
                                    f"Nitrile carbon incorporated into ring structure. Maps: {nitrile_atom_maps.intersection(ring_carbon_maps)}"
                                )
                                nitrile_to_lactone = True

                        # Check for relevant reactions even if we don't find a lactone structure yet
                        if not nitrile_to_lactone:
                            # Check reaction types that could be part of nitrile to lactone pathway
                            relevant_reactions = [
                                "Oxidation of nitrile to carboxylic acid",
                                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                "Intramolecular transesterification/Lactone formation",
                                "Esterification of Carboxylic Acids",
                            ]

                            for reaction_type in relevant_reactions:
                                if checker.check_reaction(reaction_type, rsmi):
                                    print(f"Found relevant reaction: {reaction_type}")
                                    nitrile_to_lactone = True
                                    break

                            # If we have a nitrile in reactant and an ester in product, this could be a step
                            # in the nitrile to lactone transformation
                            if not nitrile_to_lactone and nitrile_found and has_ester:
                                print(
                                    "Found nitrile to ester transformation - potential intermediate step"
                                )
                                nitrile_to_lactone = True

                            # Additional check for carboxylic acid formation from nitrile
                            if (
                                not nitrile_to_lactone
                                and nitrile_found
                                and checker.check_fg("Carboxylic acid", product)
                            ):
                                print(
                                    "Found nitrile to carboxylic acid transformation - potential intermediate step"
                                )
                                nitrile_to_lactone = True

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitrile_to_lactone
