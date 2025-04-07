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
    This function detects a convergent synthesis strategy where multiple complex fragments
    are joined in the late stages of synthesis.
    """
    has_convergent_synthesis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_synthesis

        if node["type"] == "reaction" and depth <= 3:  # Late stage (expanded depth threshold)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if we have multiple reactants
                if len(reactants) >= 2:
                    complex_reactants = 0
                    complex_reactant_smiles = []

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # Consider a reactant complex if it has many atoms OR is a key reactive fragment
                            if (
                                mol.GetNumHeavyAtoms() >= 8
                                or checker.check_fg("Formaldehyde", reactant)
                                or checker.check_fg("Aldehyde", reactant)
                            ):
                                complex_reactants += 1
                                complex_reactant_smiles.append(reactant)
                                print(
                                    f"Found complex reactant: {reactant} with {mol.GetNumHeavyAtoms()} heavy atoms"
                                )

                    # Check for coupling reaction types
                    is_coupling = False
                    coupling_reactions = [
                        "Suzuki",
                        "Negishi",
                        "Stille",
                        "Heck",
                        "Sonogashira",
                        "Buchwald-Hartwig",
                        "N-arylation",
                        "Schotten-Baumann to ester",
                        "Williamson Ether Synthesis",
                        "Acylation of Nitrogen Nucleophiles",
                        "Ugi reaction",
                        "Mitsunobu",
                        "Chan-Lam",
                        "reductive amination",
                        "Eschweiler-Clarke Primary Amine Methylation",
                        "Reductive methylation of primary amine with formaldehyde",
                    ]

                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Confirmed {rxn_type} coupling reaction")
                            is_coupling = True
                            break

                    # Check for N-methylation reactions specifically
                    is_n_methylation = False
                    if not is_coupling:
                        methylation_reactions = [
                            "N-methylation",
                            "Methylation with MeI_primary",
                            "Methylation with MeI_secondary",
                            "Methylation with DMS",
                        ]
                        for rxn_type in methylation_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Confirmed {rxn_type} reaction")
                                is_n_methylation = True
                                break

                    # Check for convergent synthesis pattern
                    # Either multiple complex reactants OR one complex reactant with a coupling/methylation reaction
                    if complex_reactants >= 2 or (
                        complex_reactants >= 1 and (is_coupling or is_n_methylation)
                    ):
                        # Verify that fragments are actually joined in the product
                        product_mol = Chem.MolFromSmiles(product)

                        if product_mol and product_mol.GetNumHeavyAtoms() > 10:
                            print(
                                f"Found convergent synthesis at depth {depth} with {complex_reactants} complex reactants"
                            )
                            print(f"Product has {product_mol.GetNumHeavyAtoms()} heavy atoms")

                            # Check if this is a simple protection/deprotection
                            simple_transformations = ["Protection", "Deprotection"]
                            is_simple = False

                            for rxn_type in simple_transformations:
                                if checker.check_reaction(rxn_type, rsmi):
                                    print(
                                        f"This is a simple {rxn_type} reaction, not convergent synthesis"
                                    )
                                    is_simple = True
                                    break

                            if not is_simple:
                                # Check if there's a significant change in the molecule
                                # For methylation reactions, check if an amine is being methylated
                                if is_n_methylation:
                                    for reactant in reactants:
                                        if checker.check_fg(
                                            "Primary amine", reactant
                                        ) or checker.check_fg("Secondary amine", reactant):
                                            print(
                                                "Confirmed amine methylation in convergent synthesis"
                                            )
                                            has_convergent_synthesis = True
                                            break
                                else:
                                    # Mark as convergent synthesis if it meets our criteria
                                    has_convergent_synthesis = True
                                    print("Confirmed convergent synthesis strategy")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {has_convergent_synthesis}")
    return has_convergent_synthesis
