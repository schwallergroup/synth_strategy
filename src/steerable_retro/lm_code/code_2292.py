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
    This function detects the construction of a multi-heterocyclic scaffold with biaryl linkages.
    """
    # The final product is the root node of the synthesis route
    final_product_smiles = route["smiles"]
    final_product = Chem.MolFromSmiles(final_product_smiles)

    if final_product is None:
        print("Could not parse final product SMILES")
        return False

    # Check for multiple heterocycles in the final product
    heterocycle_types = [
        "thiazole",
        "benzimidazole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "isoxazole",
        "pyrimidine",
        "pyrazine",
        "indole",
        "quinoline",
        "isoquinoline",
        "furan",
        "thiophene",
    ]

    heterocycles_found = []
    for ring_type in heterocycle_types:
        if checker.check_ring(ring_type, final_product_smiles):
            heterocycles_found.append(ring_type)

    print(f"Heterocycles found in final product: {heterocycles_found}")

    if len(heterocycles_found) < 2:
        print("Less than 2 heterocycles found in final product")
        return False

    # Check for biaryl linkages in the final product
    # A biaryl linkage is a bond between two aromatic rings
    ring_info = final_product.GetRingInfo()
    ring_atoms = [
        set(ring)
        for ring in ring_info.AtomRings()
        if all(final_product.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
    ]

    print(f"Found {len(ring_atoms)} aromatic rings in final product")

    # Check if there are bonds between different aromatic rings
    biaryl_bonds = []
    for bond in final_product.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()

        # Find which rings these atoms belong to
        begin_rings = [i for i, ring in enumerate(ring_atoms) if begin_atom in ring]
        end_rings = [i for i, ring in enumerate(ring_atoms) if end_atom in ring]

        # If atoms belong to different aromatic rings, it's a biaryl bond
        # Make sure the rings are truly different (not just one being a subset of the other)
        if (
            begin_rings
            and end_rings
            and not set(begin_rings).issubset(set(end_rings))
            and not set(end_rings).issubset(set(begin_rings))
        ):
            biaryl_bonds.append((begin_atom, end_atom))

    print(f"Found {len(biaryl_bonds)} biaryl bonds in final product")

    if not biaryl_bonds:
        print("No biaryl linkages found in final product")
        return False

    # Now check if the synthesis route involves heterocycle formation and biaryl coupling
    heterocycle_formation = False
    biaryl_coupling = False

    def check_reactions(node):
        nonlocal heterocycle_formation, biaryl_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]
            reactants_smiles = rxn_smiles.split(">")[0]
            product_smiles = rxn_smiles.split(">")[-1]

            # Check for heterocycle formation reactions
            # Method 1: Check if a heterocycle appears in the product but not in the reactants
            for ring in heterocycle_types:
                if not checker.check_ring(ring, reactants_smiles) and checker.check_ring(
                    ring, product_smiles
                ):
                    print(f"Found heterocycle formation for {ring}: {rxn_smiles}")
                    heterocycle_formation = True
                    break

            # Method 2: Check for named heterocycle formation reactions
            if not heterocycle_formation:
                for ring in heterocycle_types:
                    reaction_name = f"{{{ring}}}"
                    if checker.check_reaction(reaction_name, rxn_smiles):
                        print(f"Found heterocycle formation reaction for {ring}: {rxn_smiles}")
                        heterocycle_formation = True
                        break

            # Check for biaryl coupling reactions
            coupling_reactions = [
                "Suzuki",
                "{Suzuki}",
                "Stille",
                "{Stille}",
                "Negishi",
                "{Negishi}",
                "Kumada cross-coupling",
                "Hiyama-Denmark Coupling",
                "Ullmann condensation",
                "Buchwald-Hartwig",
                "{Buchwald-Hartwig}",
                "{N-arylation_heterocycles}",
                "Aryllithium cross-coupling",
                "decarboxylative_coupling",
                "{decarboxylative_coupling}",
            ]

            for rxn_type in coupling_reactions:
                if checker.check_reaction(rxn_type, rxn_smiles):
                    print(f"Found potential biaryl coupling reaction ({rxn_type}): {rxn_smiles}")

                    # Verify that the coupling creates a bond between aromatic rings
                    product = Chem.MolFromSmiles(product_smiles)
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]

                    if product and all(r is not None for r in reactants):
                        # Check if the product has a biaryl bond that doesn't exist in any reactant
                        product_ring_info = product.GetRingInfo()
                        product_ring_atoms = [
                            set(ring)
                            for ring in product_ring_info.AtomRings()
                            if all(product.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
                        ]

                        # Find biaryl bonds in product
                        product_biaryl_bonds = []
                        for bond in product.GetBonds():
                            begin_atom = bond.GetBeginAtomIdx()
                            end_atom = bond.GetEndAtomIdx()

                            begin_rings = [
                                i for i, ring in enumerate(product_ring_atoms) if begin_atom in ring
                            ]
                            end_rings = [
                                i for i, ring in enumerate(product_ring_atoms) if end_atom in ring
                            ]

                            if (
                                begin_rings
                                and end_rings
                                and not set(begin_rings).issubset(set(end_rings))
                                and not set(end_rings).issubset(set(begin_rings))
                            ):
                                print(f"Confirmed biaryl bond in product: {begin_atom}-{end_atom}")
                                biaryl_coupling = True
                                break

                    if biaryl_coupling:
                        break

        for child in node.get("children", []):
            check_reactions(child)
            if heterocycle_formation and biaryl_coupling:
                break  # Early exit if we found both

    check_reactions(route)

    # Return True if we found multiple heterocycles connected by biaryl linkages in the final product
    # AND the synthesis route involves heterocycle formation OR biaryl coupling
    result = (
        len(heterocycles_found) >= 2
        and len(biaryl_bonds) > 0
        and (heterocycle_formation or biaryl_coupling)
    )
    print(
        f"Final result: {result}, heterocycle_formation: {heterocycle_formation}, biaryl_coupling: {biaryl_coupling}"
    )
    return result
