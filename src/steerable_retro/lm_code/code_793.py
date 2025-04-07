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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects if the synthesis route involves a transition from linear to convergent approach,
    where two complex fragments are combined in a key disconnection step.
    """
    complex_fragment_combination = False

    def dfs_traverse(node):
        nonlocal complex_fragment_combination

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if we have multiple reactants combining to form a product
                if len(reactants_smiles) >= 2:
                    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol and all(r for r in reactant_mols):
                        # Count complex reactants (those with multiple rings and sufficient atoms)
                        complex_reactants = 0
                        for r in reactant_mols:
                            ring_count = r.GetRingInfo().NumRings()
                            atom_count = r.GetNumHeavyAtoms()
                            if ring_count >= 1 and atom_count >= 8:
                                complex_reactants += 1

                        # Check if this is a coupling reaction
                        is_coupling = False
                        reaction_smiles = node["metadata"].get("smiles", "")
                        coupling_reactions = [
                            "Suzuki",
                            "Negishi",
                            "Stille",
                            "Heck",
                            "Sonogashira",
                            "Buchwald-Hartwig",
                            "Ullmann",
                            "N-arylation",
                        ]

                        for rxn_type in coupling_reactions:
                            if checker.check_reaction(rxn_type, reaction_smiles):
                                is_coupling = True
                                print(f"Detected {rxn_type} coupling reaction")
                                break

                        # If we have at least 2 complex reactants combining in a coupling reaction
                        if complex_reactants >= 2 or (complex_reactants >= 1 and is_coupling):
                            complex_fragment_combination = True
                            print(
                                f"Complex fragment combination found: {complex_reactants} complex reactants"
                            )
                            print(f"Reaction SMILES: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Linear to convergent strategy: {complex_fragment_combination}")
    return complex_fragment_combination
