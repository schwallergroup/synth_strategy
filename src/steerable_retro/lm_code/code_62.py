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
    This function detects if the synthetic route employs a strategy of introducing
    heteroatoms (N, O, S, etc.) into the molecule.
    """
    # List of functional groups containing heteroatoms
    n_containing_fgs = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Amide",
        "Nitrile",
        "Nitro group",
        "Azide",
        "Aniline",
        "Pyridine",
        "Imidazole",
        "Triazole",
        "Tetrazole",
    ]

    o_containing_fgs = [
        "Alcohol",
        "Ether",
        "Ketone",
        "Aldehyde",
        "Carboxylic acid",
        "Ester",
        "Phenol",
        "Furan",
        "Epoxide",
        "Peroxide",
    ]

    s_containing_fgs = [
        "Thiol",
        "Sulfide",
        "Sulfone",
        "Sulfoxide",
        "Thiophene",
        "Thiazole",
        "Thioester",
        "Thiourea",
    ]

    other_heteroatom_fgs = [
        "Phosphate ester",
        "Boronic acid",
        "Boronic ester",
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
    ]

    # Reactions that typically introduce heteroatoms
    heteroatom_introducing_reactions = [
        "Amination",
        "Nitration",
        "Oxidation",
        "Sulfonation",
        "Halogenation",
        "Buchwald-Hartwig",
        "Chan-Lam",
        "Reductive amination",
        "Mitsunobu",
        "Williamson Ether Synthesis",
        "Azide formation",
        "Nitrile formation",
    ]

    heteroatom_introductions = 0

    def dfs_traverse(node, depth=0):
        nonlocal heteroatom_introductions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction introduces heteroatoms via functional groups
            introduced_n_fg = False
            introduced_o_fg = False
            introduced_s_fg = False
            introduced_other_fg = False

            # Check for N-containing functional groups
            for fg in n_containing_fgs:
                if not any(checker.check_fg(fg, r) for r in reactants) and checker.check_fg(
                    fg, product
                ):
                    introduced_n_fg = True
                    print(f"Found N-containing FG introduction: {fg} in reaction: {rsmi}")
                    break

            # Check for O-containing functional groups
            for fg in o_containing_fgs:
                if not any(checker.check_fg(fg, r) for r in reactants) and checker.check_fg(
                    fg, product
                ):
                    introduced_o_fg = True
                    print(f"Found O-containing FG introduction: {fg} in reaction: {rsmi}")
                    break

            # Check for S-containing functional groups
            for fg in s_containing_fgs:
                if not any(checker.check_fg(fg, r) for r in reactants) and checker.check_fg(
                    fg, product
                ):
                    introduced_s_fg = True
                    print(f"Found S-containing FG introduction: {fg} in reaction: {rsmi}")
                    break

            # Check for other heteroatom-containing functional groups
            for fg in other_heteroatom_fgs:
                if not any(checker.check_fg(fg, r) for r in reactants) and checker.check_fg(
                    fg, product
                ):
                    introduced_other_fg = True
                    print(f"Found other heteroatom FG introduction: {fg} in reaction: {rsmi}")
                    break

            # Check for specific heteroatom-introducing reactions
            reaction_introduces_heteroatom = False
            for rxn_type in heteroatom_introducing_reactions:
                try:
                    if checker.check_reaction(rxn_type, rsmi):
                        reaction_introduces_heteroatom = True
                        print(f"Found heteroatom-introducing reaction: {rxn_type} in {rsmi}")
                        break
                except:
                    continue

            # Count this as a heteroatom introduction if either condition is met
            if (
                introduced_n_fg
                or introduced_o_fg
                or introduced_s_fg
                or introduced_other_fg
                or reaction_introduces_heteroatom
            ):
                heteroatom_introductions += 1
                print(f"Total heteroatom introductions so far: {heteroatom_introductions}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if at least 2 heteroatom introductions are found
    return heteroatom_introductions >= 2
