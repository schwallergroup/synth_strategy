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
    This function detects heterocycle formation in the second half of the synthesis.
    Focuses on common heterocycles like furan, pyrrole, thiophene, etc.
    """
    heterocycle_formation_detected = False

    # List of heterocycles to check
    heterocycles = [
        "furan",
        "pyrrole",
        "thiophene",
        "oxazole",
        "thiazole",
        "imidazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
    ]

    # List of heterocycle-forming reactions to check
    heterocycle_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "pyrazole",
        "Fischer indole",
        "indole",
        "oxadiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction" and depth <= 3:  # Late stage of synthesis (lower depth)
            try:
                # Get reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a heterocycle-forming reaction
                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Heterocycle-forming reaction detected: {reaction_type}")
                        heterocycle_formation_detected = True
                        return

                # If no specific reaction type matched, check for heterocycle formation by structure
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol and all(r for r in reactants_mols):
                    # Check for heterocycle formation
                    for heterocycle in heterocycles:
                        product_has_heterocycle = checker.check_ring(heterocycle, product_part)
                        reactants_have_heterocycle = any(
                            checker.check_ring(heterocycle, r) for r in reactants_part.split(".")
                        )

                        if product_has_heterocycle and not reactants_have_heterocycle:
                            # Verify ring count change
                            reactants_ring_count = sum(
                                [r.GetRingInfo().NumRings() for r in reactants_mols]
                            )
                            product_ring_count = product_mol.GetRingInfo().NumRings()
                            print(
                                f"Ring count: Reactants={reactants_ring_count}, Product={product_ring_count}"
                            )

                            # Check if there's a ring transformation (not necessarily just an increase)
                            if product_ring_count != reactants_ring_count:
                                print(
                                    f"Heterocycle ({heterocycle}) formation detected at depth {depth}"
                                )
                                heterocycle_formation_detected = True
                                return

                            # Even if ring count doesn't change, check if we have a ring transformation
                            # This handles cases where one ring type is converted to another
                            reactant_rings = set()
                            for r_mol in reactants_mols:
                                for ring in heterocycles:
                                    if checker.check_ring(ring, Chem.MolToSmiles(r_mol)):
                                        reactant_rings.add(ring)

                            product_rings = set()
                            for ring in heterocycles:
                                if checker.check_ring(ring, Chem.MolToSmiles(product_mol)):
                                    product_rings.add(ring)

                            new_rings = product_rings - reactant_rings
                            if new_rings:
                                print(f"New heterocycle(s) formed: {new_rings}")
                                heterocycle_formation_detected = True
                                return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return heterocycle_formation_detected
