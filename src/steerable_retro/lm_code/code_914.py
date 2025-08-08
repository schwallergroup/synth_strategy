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
    This function detects if there's a ring opening/fragmentation step in the synthesis.
    In retrosynthetic analysis, ring opening appears as ring formation when moving from
    target to starting materials.
    """
    found_ring_opening = False

    # Common ring types to check for
    common_ring_types = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "thiophene",
        "benzene",
        "naphthalene",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "cycloheptane",
    ]

    # Known ring-opening reaction types
    ring_opening_reactions = [
        "Ring opening of epoxide with amine",
        "Acetal hydrolysis to diol",
        "Acetal hydrolysis to aldehyde",
        "Ketal hydrolysis to ketone",
        "Retro-Diels-Alder from oxazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_ring_opening

        # Early return if we already found a ring opening
        if found_ring_opening:
            return

        if node["type"] == "reaction":
            if node.get("metadata", {}).get("rsmi"):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a known ring-opening reaction type
                for rxn_type in ring_opening_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        found_ring_opening = True
                        print(f"Found ring opening reaction at depth {depth}: {rxn_type}")
                        print(f"Reaction SMILES: {rsmi}")
                        return

                # In retrosynthesis, the "reactants" in rsmi are actually the products of the forward reaction
                # and the "products" in rsmi are actually the reactants of the forward reaction
                reactants_part = rsmi.split(">")[0]  # Forward products (retro reactants)
                product_part = rsmi.split(">")[-1]  # Forward reactants (retro products)

                reactants = reactants_part.split(".")
                products = product_part.split(".")

                # Track rings in reactants and products
                reactant_rings = {}
                product_rings = {}

                # Process reactants (forward products)
                for reactant in reactants:
                    for ring_name in common_ring_types:
                        if checker.check_ring(ring_name, reactant):
                            if ring_name not in reactant_rings:
                                reactant_rings[ring_name] = 0
                            reactant_rings[ring_name] += 1

                # Process products (forward reactants)
                for product in products:
                    for ring_name in common_ring_types:
                        if checker.check_ring(ring_name, product):
                            if ring_name not in product_rings:
                                product_rings[ring_name] = 0
                            product_rings[ring_name] += 1

                # Check for ring opening - in retrosynthesis, this means more rings in products than reactants
                # (since products are actually the reactants in the forward direction)
                for ring_name, count in reactant_rings.items():
                    product_count = product_rings.get(ring_name, 0)
                    if count < product_count:
                        found_ring_opening = True
                        print(f"Found ring opening at depth {depth}: {ring_name} ring")
                        print(
                            f"Forward reaction: {product_count} {ring_name} rings → {count} {ring_name} rings"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        return

                # Also check total ring count as a fallback
                total_reactant_rings = sum(reactant_rings.values())
                total_product_rings = sum(product_rings.values())

                if total_reactant_rings < total_product_rings and not found_ring_opening:
                    found_ring_opening = True
                    print(f"Found general ring opening at depth {depth}")
                    print(
                        f"Forward reaction total rings: {total_product_rings} → {total_reactant_rings}"
                    )
                    print(f"Reaction SMILES: {rsmi}")
                    return

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_ring_opening
