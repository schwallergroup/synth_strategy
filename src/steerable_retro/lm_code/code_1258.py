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
    This function detects a synthetic strategy involving functionalization of an indazole core
    that is maintained throughout the synthesis.
    """
    indazole_present = False
    indazole_functionalization = False

    # Common functional groups to check
    functional_groups = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Carboxylic acid",
        "Ester",
        "Ketone",
        "Aldehyde",
        "Alkyne",
        "Nitrile",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Nitro group",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Ether",
        "Phenol",
        "Boronic acid",
        "Boronic ester",
    ]

    # Common reaction types that might functionalize indazole
    reaction_types = [
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Reduction of nitrile to amine",
        "Reduction of nitro groups to amines",
        "Esterification of Carboxylic Acids",
        "Williamson Ether Synthesis",
        "Reduction of ester to primary alcohol",
    ]

    def dfs_traverse(node):
        nonlocal indazole_present, indazole_functionalization

        if node["type"] == "mol" and node.get("smiles"):
            # Check for indazole core
            if checker.check_ring("indazole", node["smiles"]):
                indazole_present = True
                print(f"Indazole core found in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if indazole is present in both reactants and product
            product_has_indazole = checker.check_ring("indazole", product_smiles)
            reactants_with_indazole = [
                r for r in reactants_smiles if checker.check_ring("indazole", r)
            ]

            if product_has_indazole and reactants_with_indazole:
                print(f"Indazole core maintained in reaction: {rsmi}")

                # Check for known reaction types that might functionalize indazole
                reaction_type_found = False
                for rxn_type in reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected reaction type: {rxn_type}")
                        indazole_functionalization = True
                        reaction_type_found = True
                        break

                # Check for functional group changes regardless of reaction type
                product_fgs = set()
                for fg in functional_groups:
                    if checker.check_fg(fg, product_smiles):
                        product_fgs.add(fg)

                reactant_fgs = set()
                for r in reactants_with_indazole:
                    for fg in functional_groups:
                        if checker.check_fg(fg, r):
                            reactant_fgs.add(fg)

                # If there's a difference in functional groups, it's a functionalization
                if product_fgs != reactant_fgs:
                    indazole_functionalization = True
                    print(f"Indazole functionalization detected in reaction: {rsmi}")
                    print(f"Product FGs: {product_fgs}")
                    print(f"Reactant FGs: {reactant_fgs}")
                    print(f"New FGs: {product_fgs - reactant_fgs}")
                    print(f"Lost FGs: {reactant_fgs - product_fgs}")

                    # Check for specific indazole-related transformations
                    if "Ester" in product_fgs and "Primary alcohol" in reactant_fgs:
                        if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                            print("Detected esterification of indazole alcohol")

                    if "Primary alcohol" in product_fgs and "Ester" in reactant_fgs:
                        if checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                            print("Detected reduction of indazole ester to alcohol")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    strategy_detected = indazole_present and indazole_functionalization
    print(f"Indazole core present: {indazole_present}")
    print(f"Indazole functionalization detected: {indazole_functionalization}")
    print(f"Strategy detected: {strategy_detected}")

    return strategy_detected
