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
    Detects if the synthetic route involves sequential functionalization of a molecule
    containing multiple heterocyclic systems (imidazole, furan, etc.).
    """
    # List of heterocycles to check
    heterocycle_list = [
        "imidazole",
        "furan",
        "pyrrole",
        "thiophene",
        "pyridine",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "indole",
        "benzofuran",
        "benzothiophene",
        "benzimidazole",
    ]

    # List of functional groups to check for addition
    fg_list = [
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aniline",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Nitrile",
        "Nitro group",
        "Aldehyde",
        "Ketone",
        "Alcohol",
        "Phenol",
    ]

    # Track functionalized molecules and their depths
    functionalized_molecules = {}

    def dfs_traverse(node, depth=0):
        # Process molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this molecule has multiple heterocycles
            heterocycle_count = sum(
                1 for ring in heterocycle_list if checker.check_ring(ring, mol_smiles)
            )
            if heterocycle_count >= 2:
                print(f"Found molecule with {heterocycle_count} heterocycles: {mol_smiles}")

        # Process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Processing reaction: {rsmi}")

            try:
                # Find reactants with multiple heterocycles
                multi_heterocycle_reactants = []
                for r_smiles in reactants_smiles:
                    heterocycle_count = sum(
                        1 for ring in heterocycle_list if checker.check_ring(ring, r_smiles)
                    )
                    if heterocycle_count >= 2:
                        multi_heterocycle_reactants.append(r_smiles)

                if multi_heterocycle_reactants:
                    # Check if this reaction adds a new functional group
                    for reactant_smiles in multi_heterocycle_reactants:
                        for fg in fg_list:
                            # Check if the functional group is added in this reaction
                            if not checker.check_fg(fg, reactant_smiles) and checker.check_fg(
                                fg, product_smiles
                            ):
                                print(
                                    f"Found functionalization: {fg} added to multi-heterocycle system"
                                )

                                # Store the functionalized molecule with its depth
                                functionalized_molecules[product_smiles] = depth

                                # Check if this product is used as a reactant in another functionalization
                                for prev_product, prev_depth in functionalized_molecules.items():
                                    if prev_product in reactants_smiles and prev_depth != depth:
                                        print(
                                            f"Found sequential functionalization at depths {prev_depth} and {depth}"
                                        )
                                        return True
            except Exception as e:
                print(f"Error processing reaction SMILES: {rsmi}, Error: {str(e)}")

        # Recursively process children
        for child in node.get("children", []):
            result = dfs_traverse(child, depth + 1)
            if result:
                return True

        return False

    # Start traversal
    sequential_found = dfs_traverse(route)

    # If we didn't find sequential functionalization during traversal,
    # check if we have at least 2 functionalization steps at different depths
    if not sequential_found and len(functionalized_molecules) >= 2:
        depths = list(functionalized_molecules.values())
        # Check if we have at least 2 different depths
        if len(set(depths)) >= 2:
            print(
                f"Found {len(functionalized_molecules)} functionalization steps at depths {depths}"
            )
            return True

    return sequential_found
