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
    Detects if the synthesis route includes aromatization of a cyclohexadiene-like system
    to form an aromatic ring.
    """
    aromatization_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal aromatization_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            print(f"Depth {depth} - Examining reaction: {rsmi}")

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Count aromatic atoms before and after
                    reactant_aromatic_atoms = sum(
                        1 for atom in reactants_mol.GetAtoms() if atom.GetIsAromatic()
                    )
                    product_aromatic_atoms = sum(
                        1 for atom in product_mol.GetAtoms() if atom.GetIsAromatic()
                    )

                    print(
                        f"  Aromatic atoms: Reactants={reactant_aromatic_atoms}, Product={product_aromatic_atoms}"
                    )

                    # Check for aromatization (forward direction)
                    if product_aromatic_atoms > reactant_aromatic_atoms:
                        print(
                            f"  Increase in aromaticity detected: +{product_aromatic_atoms - reactant_aromatic_atoms} atoms"
                        )

                        # Check for oxidation/dehydrogenation reactions
                        is_oxidation_or_dehydrogenation = (
                            checker.check_reaction(
                                "Oxidation of aldehydes to carboxylic acids", rsmi
                            )
                            or checker.check_reaction(
                                "Oxidation of alcohol to carboxylic acid", rsmi
                            )
                            or checker.check_reaction("Dehydrogenation", rsmi)
                            or checker.check_reaction("Arene hydrogenation", rsmi)
                            or checker.check_reaction("Hydrogenation (double to single)", rsmi)
                            or checker.check_reaction("Hydrogenation (triple to double)", rsmi)
                        )

                        # Check for benzene rings in product that weren't in reactants
                        has_new_benzene = checker.check_ring(
                            "benzene", product_smiles
                        ) and not checker.check_ring("benzene", reactants_smiles)

                        if has_new_benzene:
                            print(f"  New benzene ring detected in product")

                        # Check for cyclohexadiene or cyclohexane in reactants
                        has_cyclohexadiene = checker.check_ring("cyclohexane", reactants_smiles)

                        if has_cyclohexadiene:
                            print(f"  Cyclohexane found in reactants")

                        # Confirm aromatization
                        if (has_new_benzene and has_cyclohexadiene) or (
                            is_oxidation_or_dehydrogenation and has_new_benzene
                        ):
                            print(f"  CONFIRMED: Aromatization detected in forward direction")
                            aromatization_detected = True

                        # Check for significant increase in aromaticity
                        elif product_aromatic_atoms - reactant_aromatic_atoms >= 6:
                            print(
                                f"  Significant increase in aromaticity detected: +{product_aromatic_atoms - reactant_aromatic_atoms} atoms"
                            )
                            aromatization_detected = True

                    # Check for dearomatization (reverse direction - important for retrosynthetic traversal)
                    elif reactant_aromatic_atoms > product_aromatic_atoms:
                        print(
                            f"  Decrease in aromaticity detected: -{reactant_aromatic_atoms - product_aromatic_atoms} atoms"
                        )

                        # In retrosynthetic analysis, dearomatization in forward direction = aromatization in reverse
                        # Check for benzene rings in reactants that aren't in products
                        has_benzene_in_reactants = checker.check_ring("benzene", reactants_smiles)
                        has_benzene_in_products = checker.check_ring("benzene", product_smiles)

                        if has_benzene_in_reactants and not has_benzene_in_products:
                            print(f"  Benzene ring in reactants but not in products")

                            # Check for cyclohexadiene or cyclohexane in products
                            has_cyclohexane_in_products = checker.check_ring(
                                "cyclohexane", product_smiles
                            )

                            if has_cyclohexane_in_products:
                                print(f"  Cyclohexane found in products")
                                print(
                                    f"  CONFIRMED: Dearomatization detected (aromatization in retrosynthetic direction)"
                                )
                                aromatization_detected = True

                            # Check if this is a significant dearomatization (loss of 6 or more aromatic atoms)
                            elif reactant_aromatic_atoms - product_aromatic_atoms >= 6:
                                print(
                                    f"  Significant dearomatization detected: -{reactant_aromatic_atoms - product_aromatic_atoms} atoms"
                                )
                                print(
                                    f"  CONFIRMED: Dearomatization detected (aromatization in retrosynthetic direction)"
                                )
                                aromatization_detected = True

                        # Check for specific reaction types that might involve dearomatization
                        is_dearomatization_rxn = checker.check_reaction("Arene hydrogenation", rsmi)

                        if is_dearomatization_rxn and has_benzene_in_reactants:
                            print(f"  Arene hydrogenation detected with benzene in reactants")
                            print(
                                f"  CONFIRMED: Dearomatization detected (aromatization in retrosynthetic direction)"
                            )
                            aromatization_detected = True

                    # Check for specific reaction types that might involve aromatization
                    if not aromatization_detected:
                        aromatization_rxn_types = [
                            "Dehydrogenation",
                            "Arene hydrogenation",
                            "Oxidation of alcohol to carboxylic acid",
                            "Hydrogenation (double to single)",
                            "Hydrogenation (triple to double)",
                        ]

                        for rxn_type in aromatization_rxn_types:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(
                                    f"  Potential aromatization reaction type detected: {rxn_type}"
                                )

                                # For Arene hydrogenation, we need to check in reverse direction
                                if rxn_type == "Arene hydrogenation":
                                    if checker.check_ring(
                                        "benzene", reactants_smiles
                                    ) and checker.check_ring("cyclohexane", product_smiles):
                                        print(
                                            f"  Arene hydrogenation converts benzene to cyclohexane (aromatization in reverse)"
                                        )
                                        aromatization_detected = True
                                        break
                                # For other reactions, check if benzene is involved
                                elif checker.check_ring(
                                    "benzene", product_smiles
                                ) or checker.check_ring("benzene", reactants_smiles):
                                    print(f"  Benzene ring involved in {rxn_type} reaction")
                                    aromatization_detected = True
                                    break

            except Exception as e:
                print(f"Error processing SMILES for aromatization detection: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            if not aromatization_detected:  # Stop traversal if we've already found aromatization
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: aromatization_detected = {aromatization_detected}")
    return aromatization_detected
