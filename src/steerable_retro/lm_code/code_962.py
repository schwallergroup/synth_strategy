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
    Detects nitro reduction to amine in the synthetic route.
    """
    found_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a nitro reduction reaction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Found nitro reduction reaction at depth {depth}: {rsmi}")
                found_nitro_reduction = True
                return

            # Alternative check by looking at functional group changes
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has a nitro group
            has_nitro_reactant = False
            nitro_reactant = None
            for reactant in reactants:
                if checker.check_fg("Nitro group", reactant):
                    has_nitro_reactant = True
                    nitro_reactant = reactant
                    print(f"Found reactant with nitro group: {reactant}")

                    # Check if product has an amine (any type) but no nitro group
                    has_amine_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Aniline", product)
                    )
                    has_nitro_product = checker.check_fg("Nitro group", product)

                    print(
                        f"Product has amine: {has_amine_product}, Product has nitro: {has_nitro_product}"
                    )

                    if has_amine_product and not has_nitro_product:
                        # Try to verify using atom mapping that the nitro group was converted to amine
                        try:
                            # Get atom indices for nitro group in reactant
                            nitro_indices = checker.get_fg_atom_indices(
                                "Nitro group", nitro_reactant
                            )
                            if nitro_indices:
                                # Get the atom mapping numbers for the nitro group atoms
                                reactant_mol = Chem.MolFromSmiles(nitro_reactant)
                                if reactant_mol:
                                    # Find the nitrogen atom in the nitro group
                                    for atom_indices in nitro_indices:
                                        for idx in atom_indices[
                                            0
                                        ]:  # Get the first tuple of indices
                                            atom = reactant_mol.GetAtomWithIdx(idx)
                                            if atom.GetSymbol() == "N":
                                                # Get the atom mapping number
                                                map_num = atom.GetAtomMapNum()
                                                if map_num > 0:
                                                    print(
                                                        f"Found nitro N atom with map number: {map_num}"
                                                    )
                                                    # Check if this atom is part of an amine in the product
                                                    product_mol = Chem.MolFromSmiles(product)
                                                    if product_mol:
                                                        for atom in product_mol.GetAtoms():
                                                            if (
                                                                atom.GetAtomMapNum() == map_num
                                                                and atom.GetSymbol() == "N"
                                                            ):
                                                                print(
                                                                    f"Found matching N atom in product with map number: {map_num}"
                                                                )
                                                                # Check if this N is part of an amine
                                                                found_nitro_reduction = True
                                                                print(
                                                                    f"Setting found_nitro_reduction to True"
                                                                )
                                                                return
                        except Exception as e:
                            print(f"Error in atom mapping check: {e}")

                        # If atom mapping check failed but we still have evidence of nitro reduction
                        print(f"Found nitro reduction at depth {depth}: {rsmi}")
                        print(f"Setting found_nitro_reduction to True")
                        found_nitro_reduction = True
                        return

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_nitro_reduction}")
    return found_nitro_reduction
