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
    Checks if the route contains Sonogashira couplings that build a diarylacetylene scaffold.
    """
    # Track if we've found a diarylacetylene scaffold built by Sonogashira
    found_diarylacetylene = False
    diarylacetylene_product = None

    def dfs(node, depth=0):
        nonlocal found_diarylacetylene, diarylacetylene_product

        # For reaction nodes, check if it's a Sonogashira coupling
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check if this is a Sonogashira coupling
            is_sonogashira = (
                checker.check_reaction("Sonogashira acetylene_aryl halide", rxn_smiles)
                or checker.check_reaction("Sonogashira alkyne_aryl halide", rxn_smiles)
                or checker.check_reaction("Sonogashira acetylene_aryl OTf", rxn_smiles)
                or checker.check_reaction("Sonogashira alkyne_aryl OTf", rxn_smiles)
            )

            if is_sonogashira:
                print(f"Found Sonogashira coupling at depth {depth}: {rxn_smiles}")

                # Check if the product contains a diarylacetylene scaffold
                try:
                    product = rxn_smiles.split(">")[-1]

                    # Check for alkyne group
                    has_alkyne = checker.check_fg("Alkyne", product)

                    # Check if the alkyne connects two aromatic rings
                    mol = Chem.MolFromSmiles(product)
                    if mol and has_alkyne:
                        # Find alkyne bonds
                        alkyne_bonds = [
                            bond.GetIdx()
                            for bond in mol.GetBonds()
                            if bond.GetBondType() == Chem.BondType.TRIPLE
                        ]

                        for bond_idx in alkyne_bonds:
                            bond = mol.GetBondWithIdx(bond_idx)
                            begin_atom = bond.GetBeginAtomIdx()
                            end_atom = bond.GetEndAtomIdx()

                            # Check if both atoms connected to the alkyne are part of aromatic rings
                            begin_aromatic = False
                            end_aromatic = False

                            for ring in mol.GetSSSR():
                                ring_atoms = set(ring)
                                if (
                                    begin_atom in ring_atoms
                                    and mol.GetAtomWithIdx(begin_atom).GetIsAromatic()
                                ):
                                    begin_aromatic = True
                                if (
                                    end_atom in ring_atoms
                                    and mol.GetAtomWithIdx(end_atom).GetIsAromatic()
                                ):
                                    end_aromatic = True

                            if begin_aromatic and end_aromatic:
                                found_diarylacetylene = True
                                diarylacetylene_product = product
                                print(f"Confirmed diarylacetylene scaffold in product: {product}")
                                break
                except Exception as e:
                    print(f"Error analyzing Sonogashira product: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal from the root
    dfs(route)
    return found_diarylacetylene, diarylacetylene_product
