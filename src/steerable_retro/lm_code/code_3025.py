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
    Detects if the synthesis uses an early-stage ring formation strategy.

    Early-stage means the ring formation occurs at the beginning of the synthesis,
    which corresponds to higher depth values in the retrosynthetic tree.
    """
    early_ring_formation = False
    final_product_smiles = route["smiles"] if route["type"] == "mol" else ""
    formed_rings_in_synthesis = set()

    # List of common ring formation reactions
    ring_formation_reactions = [
        "Diels-Alder",
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "Pictet-Spengler",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "pyrazole",
    ]

    # List of common ring types to check
    ring_types = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "furan",
        "thiophene",
        "benzene",
        "naphthalene",
        "indole",
        "quinoline",
        "isoquinoline",
        "pyrimidine",
        "pyrazine",
        "triazole",
        "tetrazole",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "piperidine",
        "morpholine",
        "pyrrolidine",
    ]

    # Get the final product rings
    final_product_rings = set()
    if final_product_smiles:
        final_product_mol = Chem.MolFromSmiles(final_product_smiles)
        if final_product_mol:
            for ring_type in ring_types:
                if checker.check_ring(ring_type, final_product_smiles):
                    final_product_rings.add(ring_type)

    print(f"Final product rings: {final_product_rings}")

    def dfs_traverse(node, depth=0):
        nonlocal early_ring_formation, formed_rings_in_synthesis

        if node["type"] == "reaction" and depth > 1:  # Early-stage in retrosynthetic analysis
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if valid molecules can be created
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                reactant_mols = [mol for mol in reactant_mols if mol is not None]
                product_mol = Chem.MolFromSmiles(product)

                if not product_mol or not reactant_mols:
                    return

                # Count rings in reactants and product using RingInfo
                reactant_rings = sum(mol.GetRingInfo().NumRings() for mol in reactant_mols)
                product_rings = product_mol.GetRingInfo().NumRings()

                print(f"Rings in reactants: {reactant_rings}, Rings in product: {product_rings}")

                # Check for ring formation
                ring_formation_detected = False

                # 1. Check if this is a known ring formation reaction
                is_ring_formation_rxn = False
                for rxn_type in ring_formation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_ring_formation_rxn = True
                        print(f"Ring formation reaction detected: {rxn_type}")
                        ring_formation_detected = True
                        break

                # 2. Check if product has more rings than reactants combined
                if product_rings > reactant_rings:
                    ring_formation_detected = True
                    print("Ring formation detected: product has more rings than reactants")

                # 3. Check which specific rings are formed
                formed_rings = []
                for ring_type in ring_types:
                    # Check if ring exists in product but not in any reactant
                    if checker.check_ring(ring_type, product):
                        if not any(
                            checker.check_ring(ring_type, r)
                            for r in reactants
                            if Chem.MolFromSmiles(r)
                        ):
                            formed_rings.append(ring_type)
                            formed_rings_in_synthesis.add(ring_type)
                            print(f"New ring formed: {ring_type}")

                if formed_rings:
                    ring_formation_detected = True

                # 4. Check for ring transformations (same count but different types)
                if (
                    not ring_formation_detected
                    and product_rings == reactant_rings
                    and product_rings > 0
                ):
                    product_ring_types = set()
                    reactant_ring_types = set()

                    for ring_type in ring_types:
                        if checker.check_ring(ring_type, product):
                            product_ring_types.add(ring_type)

                        for r in reactants:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol and checker.check_ring(ring_type, r):
                                reactant_ring_types.add(ring_type)

                    # If there are new ring types in the product
                    new_ring_types = product_ring_types - reactant_ring_types
                    if new_ring_types:
                        formed_rings.extend(list(new_ring_types))
                        formed_rings_in_synthesis.update(new_ring_types)
                        ring_formation_detected = True
                        print(f"Ring transformation detected: {new_ring_types}")

                # 5. Check for cyclization reactions that might not be in our predefined list
                if not ring_formation_detected and product_rings > reactant_rings:
                    # This is a catch-all for any cyclization not caught by previous checks
                    ring_formation_detected = True
                    print("Generic cyclization detected")

                if ring_formation_detected:
                    print(f"Found early-stage ring formation at depth {depth}")
                    print(f"Reaction: {rsmi}")
                    print(
                        f"Rings in reactants: {reactant_rings}, Rings in product: {product_rings}"
                    )
                    if formed_rings:
                        print(f"Rings formed: {', '.join(formed_rings)}")

                    # Set early_ring_formation to True if we detect ring formation at an early stage
                    # We don't need to check if the formed rings persist in the final product
                    early_ring_formation = True

                    # For debugging, show if formed rings are in final product
                    if final_product_rings and formed_rings_in_synthesis:
                        common_rings = final_product_rings.intersection(formed_rings_in_synthesis)
                        if common_rings:
                            print(
                                f"Formed rings that persist in final product: {', '.join(common_rings)}"
                            )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Early ring formation detected: {early_ring_formation}")
    print(f"Rings formed during synthesis: {formed_rings_in_synthesis}")

    return early_ring_formation
