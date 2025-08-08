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
    This function detects a synthetic strategy involving halogenation
    (chlorination, bromination, fluorination, or iodination).
    """
    has_halogenation = False

    # List of halogenation reaction types to check
    halogenation_reactions = [
        # Chlorination reactions
        "Aromatic chlorination",
        "Chlorination",
        "Alcohol to chloride_SOCl2",
        "Alcohol to chloride_CHCl3",
        "Alcohol to chloride_CH2Cl2",
        "Alcohol to chloride_PCl5_ortho",
        "Alcohol to chloride_POCl3_ortho",
        "Alcohol to chloride_POCl3_para",
        "Alcohol to chloride_POCl3",
        "Alcohol to chloride_HCl",
        "Alcohol to chloride_Salt",
        "Alcohol to chloride_Other",
        "Alcohol to chloride_sulfonyl chloride",
        # Bromination reactions
        "Aromatic bromination",
        "Bromination",
        "Wohl-Ziegler bromination benzyl primary",
        "Wohl-Ziegler bromination benzyl secondary",
        "Wohl-Ziegler bromination benzyl tertiary",
        "Wohl-Ziegler bromination allyl primary",
        "Wohl-Ziegler bromination allyl secondary",
        "Wohl-Ziegler bromination allyl tertiary",
        "Wohl-Ziegler bromination carbonyl primary",
        "Wohl-Ziegler bromination carbonyl secondary",
        "Wohl-Ziegler bromination carbonyl tertiary",
        # Fluorination reactions
        "Aromatic fluorination",
        "Fluorination",
        # Iodination reactions
        "Aromatic iodination",
        "Iodination",
    ]

    # Halogen functional groups to check
    halogen_fgs = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Alkenyl halide",
        "Haloalkyne",
    ]

    def contains_halogen_atom(smiles):
        """Check if a molecule contains any halogen atoms (F, Cl, Br, I)"""
        if not smiles:
            return False
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                for atom in mol.GetAtoms():
                    atomic_num = atom.GetAtomicNum()
                    # F(9), Cl(17), Br(35), I(53)
                    if atomic_num in [9, 17, 35, 53]:
                        return True
            return False
        except:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal has_halogenation

        print(f"Traversing node at depth {depth}: {node.get('smiles', 'reaction node')}")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")
                products = products_part.split(".")

                print(f"Checking reaction: {rsmi}")

                # Check if this is a known halogenation reaction
                for rxn_type in halogenation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found halogenation reaction: {rxn_type}")
                        has_halogenation = True

                # Check for halogen atoms in reactants and products
                reactants_have_halogen = any(contains_halogen_atom(r) for r in reactants)
                products_have_halogen = any(contains_halogen_atom(p) for p in products)

                # Check for halogen functional groups in reactants and products
                reactants_have_halogen_fg = any(
                    checker.check_fg(fg, reactant) for fg in halogen_fgs for reactant in reactants
                )

                products_have_halogen_fg = any(
                    checker.check_fg(fg, product) for fg in halogen_fgs for product in products
                )

                # If halogen is introduced (not present in reactants but present in products)
                if (not reactants_have_halogen and products_have_halogen) or (
                    not reactants_have_halogen_fg and products_have_halogen_fg
                ):
                    print(f"Found halogenation: halogen introduced in reaction")
                    has_halogenation = True

                # If halogen is maintained (present in both reactants and products)
                if (reactants_have_halogen and products_have_halogen) or (
                    reactants_have_halogen_fg and products_have_halogen_fg
                ):
                    print(f"Halogen maintained through reaction")
                    # This indicates halogen is part of the strategy
                    has_halogenation = True

        elif node["type"] == "mol":
            # Check if the molecule itself contains halogen atoms
            if contains_halogen_atom(node["smiles"]):
                print(f"Found molecule with halogen atom: {node['smiles']}")
                has_halogenation = True

            # Check if the molecule contains halogen functional groups
            for fg in halogen_fgs:
                if checker.check_fg(fg, node["smiles"]):
                    print(f"Found molecule with {fg}: {node['smiles']}")
                    has_halogenation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Has halogenation strategy: {has_halogenation}")

    return has_halogenation
