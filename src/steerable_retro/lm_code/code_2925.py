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
    This function detects lactam ring formation in the synthetic route.
    """
    lactam_formed = False

    def dfs_traverse(node):
        nonlocal lactam_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Split reactants into individual molecules
                    reactants_list = reactants_smiles.split(".")

                    # Check for lactam rings in product
                    lactam_in_product = False
                    if checker.check_ring("pyrrolidone", product_smiles):
                        lactam_in_product = True
                        print(f"Pyrrolidone (5-membered lactam) found in product: {product_smiles}")

                    # Check for other common lactam-containing rings
                    if not lactam_in_product and checker.check_fg(
                        "Secondary amide", product_smiles
                    ):
                        # Check if the amide is part of a ring
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            # Check if the product has a ring that contains both N and C=O
                            ring_info = product_mol.GetRingInfo()
                            if ring_info.NumRings() > 0:
                                for atom in product_mol.GetAtoms():
                                    if atom.GetSymbol() == "N" and atom.IsInRing():
                                        # Check if this nitrogen is part of an amide in a ring
                                        for neighbor in atom.GetNeighbors():
                                            if neighbor.GetSymbol() == "C":
                                                for bond in neighbor.GetBonds():
                                                    other_atom = bond.GetOtherAtom(neighbor)
                                                    if (
                                                        other_atom.GetSymbol() == "O"
                                                        and bond.GetBondType()
                                                        == Chem.BondType.DOUBLE
                                                    ):
                                                        if ring_info.AreAtomsInSameRing(
                                                            atom.GetIdx(), neighbor.GetIdx()
                                                        ):
                                                            lactam_in_product = True
                                                            print(
                                                                f"Lactam structure found in product: {product_smiles}"
                                                            )
                                                            break

                    # If lactam is in product, check if it's in any reactant
                    if lactam_in_product:
                        lactam_in_reactants = False
                        for reactant in reactants_list:
                            if checker.check_ring("pyrrolidone", reactant):
                                lactam_in_reactants = True
                                print(f"Pyrrolidone found in reactant: {reactant}")
                                break

                            # Check for other lactam structures in reactants
                            if checker.check_fg("Secondary amide", reactant):
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    ring_info = reactant_mol.GetRingInfo()
                                    if ring_info.NumRings() > 0:
                                        for atom in reactant_mol.GetAtoms():
                                            if atom.GetSymbol() == "N" and atom.IsInRing():
                                                for neighbor in atom.GetNeighbors():
                                                    if neighbor.GetSymbol() == "C":
                                                        for bond in neighbor.GetBonds():
                                                            other_atom = bond.GetOtherAtom(neighbor)
                                                            if (
                                                                other_atom.GetSymbol() == "O"
                                                                and bond.GetBondType()
                                                                == Chem.BondType.DOUBLE
                                                            ):
                                                                if ring_info.AreAtomsInSameRing(
                                                                    atom.GetIdx(), neighbor.GetIdx()
                                                                ):
                                                                    lactam_in_reactants = True
                                                                    print(
                                                                        f"Lactam structure found in reactant: {reactant}"
                                                                    )
                                                                    break

                        # If lactam is in product but not in reactants, it was formed in this reaction
                        if not lactam_in_reactants:
                            lactam_formed = True
                            print(f"Lactam formation detected in reaction: {rsmi}")

                        # Also check specific reaction types that might form lactams
                        if (
                            checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                            )
                            or checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                rsmi,
                            )
                            or checker.check_reaction(
                                "Intramolecular transesterification/Lactone formation", rsmi
                            )
                        ):
                            if not lactam_in_reactants and lactam_in_product:
                                lactam_formed = True
                                print(
                                    f"Lactam formation via specific reaction type detected: {rsmi}"
                                )

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return lactam_formed
