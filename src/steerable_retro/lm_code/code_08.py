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
    This function detects a synthetic strategy involving the reduction of an aromatic
    nitro group to an amine.
    """
    found_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if this is a nitro reduction reaction
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print(f"Found nitro reduction reaction: {rsmi}")
                        # Verify nitro group is on aromatic ring
                        for reactant in reactants:
                            if checker.check_fg("Nitro group", reactant):
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    nitro_indices = (
                                        checker.get_fg_atom_indices("Nitro group", reactant) or []
                                    )
                                    for indices_group in nitro_indices:
                                        # Get nitrogen atom of nitro group
                                        n_idx = indices_group[0]
                                        for neighbor in reactant_mol.GetAtomWithIdx(
                                            n_idx
                                        ).GetNeighbors():
                                            if (
                                                neighbor.GetIsAromatic()
                                                and neighbor.GetSymbol() == "C"
                                            ):
                                                found_nitro_reduction = True
                                                print(f"Confirmed aromatic nitro reduction: {rsmi}")
                                                break
                    else:
                        # Alternative check for nitro reduction
                        for reactant in reactants:
                            if checker.check_fg("Nitro group", reactant) and checker.check_fg(
                                "Primary amine", product
                            ):
                                # Verify it's attached to an aromatic ring
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    # Get indices of nitro groups
                                    nitro_indices = (
                                        checker.get_fg_atom_indices("Nitro group", reactant) or []
                                    )
                                    for indices_group in nitro_indices:
                                        # Get nitrogen atom of nitro group
                                        n_idx = indices_group[0]
                                        for neighbor in reactant_mol.GetAtomWithIdx(
                                            n_idx
                                        ).GetNeighbors():
                                            if (
                                                neighbor.GetIsAromatic()
                                                and neighbor.GetSymbol() == "C"
                                            ):
                                                # Check if product has amine
                                                product_mol = Chem.MolFromSmiles(product)
                                                if product_mol:
                                                    print(
                                                        f"Found potential aromatic nitro reduction: {rsmi}"
                                                    )
                                                    found_nitro_reduction = True
                                                    break
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_reduction
