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
    Detects if the synthesis includes an epoxide ring opening step.
    """
    has_ring_opening = False

    def dfs_traverse(node):
        nonlocal has_ring_opening

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Split into individual molecules
            reactant_mols = reactants_part.split(".")
            product_mols = products_part.split(".")

            # Check for epoxide in reactants
            for reactant in reactant_mols:
                if checker.check_ring("oxirane", reactant):
                    print(f"Found epoxide in reactant: {reactant}")

                    # Check if this is a known epoxide ring opening reaction
                    epoxide_reactions = [
                        "Ring opening of epoxide with amine",
                        "Oxirane functionalization with ketones",
                        "Oxirane functionalization with alkyl iodide",
                        "Oxirane functionalization with silyl chloride",
                    ]

                    for rxn_type in epoxide_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Found epoxide ring opening: {rxn_type}")
                            has_ring_opening = True
                            return

                    # If specific reaction check fails, check if epoxide is absent in products
                    epoxide_found_in_products = False
                    for product in product_mols:
                        if checker.check_ring("oxirane", product):
                            epoxide_found_in_products = True
                            break

                    if not epoxide_found_in_products:
                        # Check for C-O bonds in products (indicating ring opening)
                        for product in product_mols:
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                # Look for C-O single bonds
                                for bond in product_mol.GetBonds():
                                    begin_atom = bond.GetBeginAtom()
                                    end_atom = bond.GetEndAtom()
                                    if bond.GetBondType() == Chem.BondType.SINGLE and (
                                        (
                                            begin_atom.GetSymbol() == "C"
                                            and end_atom.GetSymbol() == "O"
                                        )
                                        or (
                                            begin_atom.GetSymbol() == "O"
                                            and end_atom.GetSymbol() == "C"
                                        )
                                    ):
                                        print(
                                            f"Found epoxide ring opening (general case): C-O bond in product"
                                        )
                                        has_ring_opening = True
                                        return

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_ring_opening
