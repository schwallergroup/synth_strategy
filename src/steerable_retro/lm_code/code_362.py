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
    This function detects a synthetic strategy involving heterocycle transformation,
    specifically heterocycle ring opening, formation, or significant modification.
    """
    heterocycle_transformed = False

    # List of heterocycles to check
    heterocycles = [
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "isothiazole",
        "pyrrole",
        "imidazole",
        "pyridine",
        "furan",
        "thiophene",
    ]

    # List of functional groups that might indicate heterocycle modification
    modification_fgs = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Nitrile",
        "Carboxylic acid",
        "Ester",
        "Aldehyde",
        "Ketone",
    ]

    # Reactions commonly used for heterocycle transformations
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
        "indole",
        "oxadiazole",
        "Fischer indole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_transformed

        if node["type"] == "reaction" and not heterocycle_transformed:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a known heterocycle-forming/transforming reaction
            for reaction_type in heterocycle_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected heterocycle transformation reaction: {reaction_type}")
                    heterocycle_transformed = True
                    return

            # Check for heterocycles in reactants and their presence/absence in product
            for heterocycle in heterocycles:
                # Check for heterocycle disappearance (transformation/opening)
                heterocycle_in_reactants = False
                heterocycle_reactant = None

                for reactant in reactants_smiles:
                    if checker.check_ring(heterocycle, reactant):
                        heterocycle_in_reactants = True
                        heterocycle_reactant = reactant
                        print(f"Found {heterocycle} in reactant: {reactant}")
                        break

                # Check for heterocycle appearance (formation)
                heterocycle_in_product = checker.check_ring(heterocycle, product_smiles)

                # Case 1: Heterocycle disappearance
                if heterocycle_in_reactants and not heterocycle_in_product:
                    heterocycle_transformed = True
                    print(
                        f"{heterocycle} ring transformation detected: present in reactant but absent in product"
                    )
                    return

                # Case 2: Heterocycle formation
                if not heterocycle_in_reactants and heterocycle_in_product:
                    heterocycle_transformed = True
                    print(
                        f"{heterocycle} ring formation detected: absent in reactants but present in product"
                    )
                    return

                # Case 3: Heterocycle modification (same type but modified structure)
                if heterocycle_in_reactants and heterocycle_in_product:
                    reactant_mol = Chem.MolFromSmiles(heterocycle_reactant)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactant_mol and product_mol:
                        # Get atom indices for the heterocycle in both molecules
                        reactant_indices = checker.get_ring_atom_indices(
                            heterocycle, heterocycle_reactant
                        )
                        product_indices = checker.get_ring_atom_indices(heterocycle, product_smiles)

                        if reactant_indices and product_indices:
                            # If the number of matches differs, consider it a transformation
                            if len(reactant_indices) != len(product_indices):
                                heterocycle_transformed = True
                                print(
                                    f"{heterocycle} ring transformation detected: different number of matches"
                                )
                                return

                            # Check for functional group changes on the heterocycle
                            for fg in modification_fgs:
                                fg_in_reactant = checker.check_fg(fg, heterocycle_reactant)
                                fg_in_product = checker.check_fg(fg, product_smiles)

                                if fg_in_reactant != fg_in_product:
                                    heterocycle_transformed = True
                                    print(
                                        f"{heterocycle} ring modification detected: change in {fg} functional group"
                                    )
                                    return

        # Traverse children (prioritize late-stage transformations by checking lower depths first)
        for child in sorted(
            node.get("children", []), key=lambda x: 0 if x["type"] == "reaction" else 1
        ):
            if not heterocycle_transformed:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Heterocycle transformation strategy detected: {heterocycle_transformed}")
    return heterocycle_transformed
