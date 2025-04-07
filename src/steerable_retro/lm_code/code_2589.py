#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a strategy involving C-C bond disconnection at a benzylic position.
    """
    has_cc_disconnection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_cc_disconnection

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1].split(".")

                # In retrosynthetic analysis, we're looking at reactions where
                # one product is formed from multiple reactants (C-C bond formation)
                if len(reactants) >= 2 and len(product) == 1:
                    product_mol = Chem.MolFromSmiles(product[0])

                    # Check if the product has a benzylic carbon
                    # This includes various benzylic patterns
                    benzylic_patterns = [
                        "[c;R]-[CH2;!R]",  # Benzylic CH2
                        "[c;R]-[CH;!R]",  # Benzylic CH
                        "[c;R]-[C;!R]",  # Benzylic C
                    ]

                    for pattern in benzylic_patterns:
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts(pattern)
                        ):
                            print(f"Found potential benzylic position in product: {product[0]}")

                            # Check if reactants contain carbonyl compounds (aldehydes, ketones)
                            has_carbonyl = False
                            for reactant in reactants:
                                if checker.check_fg("Aldehyde", reactant) or checker.check_fg(
                                    "Ketone", reactant
                                ):
                                    has_carbonyl = True
                                    print(f"Found carbonyl group in reactant: {reactant}")

                            # Check if this is a known C-C bond forming reaction
                            cc_forming_reactions = [
                                "Aldol condensation",
                                "Grignard_carbonyl",
                                "Wittig",
                                "Suzuki",
                                "Negishi",
                                "Heck_terminal_vinyl",
                                "Stille",
                            ]

                            for rxn_type in cc_forming_reactions:
                                if checker.check_reaction(rxn_type, rsmi):
                                    print(f"Detected {rxn_type} reaction: {rsmi}")
                                    has_cc_disconnection = True
                                    return

                            # If no specific reaction detected, check for general C-C bond formation
                            # by examining if an aromatic ring is present in both product and one reactant
                            if has_carbonyl:
                                aromatic_in_product = False
                                aromatic_in_reactant = False

                                if checker.check_ring("benzene", product[0]) or checker.check_ring(
                                    "thiophene", product[0]
                                ):
                                    aromatic_in_product = True

                                for reactant in reactants:
                                    if checker.check_ring(
                                        "benzene", reactant
                                    ) or checker.check_ring("thiophene", reactant):
                                        aromatic_in_reactant = True

                                if aromatic_in_product and aromatic_in_reactant:
                                    print(
                                        f"Detected C-C bond disconnection at benzylic position: {rsmi}"
                                    )
                                    has_cc_disconnection = True
                                    return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_cc_disconnection
