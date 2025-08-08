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
    Detects a linear synthesis strategy where a halogen substituent
    (typically on an aromatic ring) is preserved throughout the synthesis.
    """
    # Track if synthesis is linear and if halogen is preserved
    is_linear = True
    halogen_present = False

    # Track halogen atoms by their atom mapping
    halogen_atom_maps = set()

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, halogen_present

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for halogen on aromatic ring in molecules
            if checker.check_fg("Aromatic halide", mol_smiles):
                halogen_present = True
                print(f"Depth {depth}: Found aromatic halide in molecule: {mol_smiles}")

                # Track atom mappings of halogens in starting materials
                if node.get("in_stock", False):
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        for atom in mol.GetAtoms():
                            if atom.GetSymbol() in ["F", "Cl", "Br", "I"] and atom.GetIsAromatic():
                                if atom.HasProp("molAtomMapNumber"):
                                    map_num = atom.GetProp("molAtomMapNumber")
                                    halogen_atom_maps.add(map_num)
                                    print(
                                        f"Depth {depth}: Found mapped aromatic halide atom {atom.GetSymbol()} with map number {map_num} in starting material"
                                    )

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            agents_part = rsmi.split(">")[1]
            product_part = rsmi.split(">")[2]

            reactants = reactants_part.split(".")
            product = product_part

            # Check if this is a convergent step (more than one significant reactant)
            significant_reactants = [r for r in reactants if not is_likely_reagent(r)]

            # Only consider it non-linear if multiple significant reactants AND
            # more than one contains an aromatic halide
            halogen_containing_reactants = [
                r for r in significant_reactants if checker.check_fg("Aromatic halide", r)
            ]
            if len(significant_reactants) > 1 and len(halogen_containing_reactants) > 1:
                is_linear = False
                print(
                    f"Depth {depth}: Non-linear step detected with {len(halogen_containing_reactants)} halogen-containing reactants"
                )

            # Check for halogen preservation in this reaction
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Track halogen atom mappings in product
                product_has_halogen = checker.check_fg("Aromatic halide", product)
                if product_has_halogen:
                    for atom in product_mol.GetAtoms():
                        if atom.GetSymbol() in ["F", "Cl", "Br", "I"] and atom.GetIsAromatic():
                            if atom.HasProp("molAtomMapNumber"):
                                map_num = atom.GetProp("molAtomMapNumber")
                                halogen_atom_maps.add(map_num)
                                print(
                                    f"Depth {depth}: Found mapped aromatic halide atom {atom.GetSymbol()} with map number {map_num} in product"
                                )

                # Track halogen atom mappings in reactants
                for reactant in significant_reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and checker.check_fg("Aromatic halide", reactant):
                        for atom in reactant_mol.GetAtoms():
                            if atom.GetSymbol() in ["F", "Cl", "Br", "I"] and atom.GetIsAromatic():
                                if atom.HasProp("molAtomMapNumber"):
                                    map_num = atom.GetProp("molAtomMapNumber")
                                    halogen_atom_maps.add(map_num)
                                    print(
                                        f"Depth {depth}: Found mapped aromatic halide atom {atom.GetSymbol()} with map number {map_num} in reactant"
                                    )

        # Process children (reactants in retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

        return True

    def is_likely_reagent(smiles):
        """Helper function to identify common reagents that shouldn't count for convergent synthesis"""
        # Common reagents and small molecules
        common_reagents = [
            "O",
            "N",
            "C",
            "CC",
            "[NH3]",
            "CO",
            "CCO",
            "CN",
            "H2O",
            "HCl",
            "NaOH",
            "KOH",
            "H2SO4",
        ]

        # Check if it's in our list of common reagents
        if smiles in common_reagents:
            return True

        # Check if it's a very small molecule (likely a reagent) but not if it contains an aromatic halide
        if len(smiles) < 4 and not checker.check_fg("Aromatic halide", smiles):
            return True

        # Try to identify common reagent patterns
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Very small molecules with few atoms are likely reagents, unless they contain an aromatic halide
            if mol.GetNumAtoms() < 3 and not checker.check_fg("Aromatic halide", smiles):
                return True

            # Check for common acid/base reagents
            if (
                checker.check_fg("Carboxylic acid", smiles)
                and mol.GetNumAtoms() < 5
                and not checker.check_fg("Aromatic halide", smiles)
            ):
                return True

        return False

    # Start traversal
    dfs_traverse(route)

    # Check if we found any halogen atoms and if the synthesis is linear
    result = is_linear and halogen_present
    print(f"Final result: is_linear={is_linear}, halogen_present={halogen_present}")
    return result
