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
    Detects if the synthesis route involves a Suzuki coupling (C-C bond formation
    between an aryl/heteroaryl halide and a boronic acid/ester).
    """
    has_suzuki = False

    def dfs_traverse(node):
        nonlocal has_suzuki

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            boronic_pattern = re.compile(r"B\(O\)|OB\(O\)")
            has_boronic = any(boronic_pattern.search(r) for r in reactants_smiles)

            # Check for aryl/heteroaryl halide pattern
            halide_pattern = re.compile(r"Br|I|Cl")
            has_halide = any(halide_pattern.search(r) for r in reactants_smiles)

            if has_boronic and has_halide:
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(r for r in reactants):
                        # Check if a new C-C bond is formed between aromatic carbons
                        aromatic_c_c_bonds_product = set()
                        for bond in product.GetBonds():
                            if (
                                bond.GetBeginAtom().GetAtomicNum() == 6
                                and bond.GetEndAtom().GetAtomicNum() == 6
                                and bond.GetBeginAtom().GetIsAromatic()
                                and bond.GetEndAtom().GetIsAromatic()
                            ):
                                idx1, idx2 = sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
                                aromatic_c_c_bonds_product.add((idx1, idx2))

                        aromatic_c_c_bonds_reactants = set()
                        for r in reactants:
                            for bond in r.GetBonds():
                                if (
                                    bond.GetBeginAtom().GetAtomicNum() == 6
                                    and bond.GetEndAtom().GetAtomicNum() == 6
                                    and bond.GetBeginAtom().GetIsAromatic()
                                    and bond.GetEndAtom().GetIsAromatic()
                                ):
                                    idx1, idx2 = sorted(
                                        [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
                                    )
                                    aromatic_c_c_bonds_reactants.add((idx1, idx2))

                        if len(aromatic_c_c_bonds_product) > len(aromatic_c_c_bonds_reactants):
                            print("Suzuki coupling detected")
                            has_suzuki = True
                except:
                    print("Error processing reaction SMILES for Suzuki detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_suzuki
