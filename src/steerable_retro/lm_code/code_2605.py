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
    Detects if the synthesis route involves formation of a heterocycle
    (such as thiazole, oxazole, imidazole, etc.) in the final step of the synthesis.
    """
    heterocycle_formed = False

    # List of heterocycles to check
    heterocycles = [
        "thiazole",
        "oxazole",
        "isoxazole",
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "furan",
        "pyrrole",
        "benzothiazole",
        "benzoxazole",
        "benzimidazole",
        "oxadiazole",
        "thiadiazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
    ]

    # Corresponding reaction types
    reaction_types = [
        "{thiazole}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "{pyrazole}",
        "{oxadiazole}",
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "{Paal-Knorr pyrrole}",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "{Fischer indole}",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "{Huisgen_Cu-catalyzed_1,4-subst}",
        "{Huisgen_Ru-catalyzed_1,5_subst}",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
    ]

    print(f"Starting analysis for late-stage heterocycle formation")

    def check_heterocycle_formation(node, depth=0):
        nonlocal heterocycle_formed

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate reaction (late-stage)
            try:
                # Get reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing late-stage reaction at depth {depth}: {rsmi}")

                # Check for heterocycle formation
                for heterocycle in heterocycles:
                    # Check if heterocycle is in product
                    if checker.check_ring(heterocycle, product_smiles):
                        print(f"Found {heterocycle} in product")

                        # Check if heterocycle is not in any reactant
                        heterocycle_in_reactants = any(
                            checker.check_ring(heterocycle, r_smiles)
                            for r_smiles in reactants_smiles
                        )

                        if not heterocycle_in_reactants:
                            print(f"{heterocycle} not found in reactants, checking reaction type")

                            # Check for specific reaction types
                            for reaction_type in reaction_types:
                                if checker.check_reaction(reaction_type, rsmi):
                                    heterocycle_formed = True
                                    print(
                                        f"Detected {reaction_type} formation in late-stage step (depth {depth})"
                                    )
                                    return  # Exit once we find a match

                            # If no specific reaction type matched but heterocycle is formed
                            print(
                                f"No specific reaction type matched, but heterocycle appears to be formed"
                            )

                            # Since we've verified a heterocycle appears in product but not reactants,
                            # it's likely a heterocycle formation even if we don't match a specific reaction type
                            heterocycle_formed = True
                            return
                        else:
                            print(
                                f"{heterocycle} also found in reactants, not a formation reaction"
                            )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            if not heterocycle_formed:  # Stop traversal if we've already found a match
                check_heterocycle_formation(child, depth + 1)

    # Start traversal
    check_heterocycle_formation(route)

    print(f"Heterocycle formation detected: {heterocycle_formed}")
    return heterocycle_formed
