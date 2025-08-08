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
    This function detects a strategy where nitro reduction is performed
    to create an amine that is later used in a cyclization reaction.
    """
    nitro_reduction_found = False
    subsequent_cyclization_found = False
    reaction_depths = {}
    nitro_reduction_depth = None
    reduced_nitro_mol = None

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, subsequent_cyclization_found, nitro_reduction_depth, reduced_nitro_mol

        if node["type"] == "reaction":
            # Store reaction at its depth
            reaction_depths[depth] = node

            # Check if this is a nitro reduction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction pattern
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Nitro reduction detected at depth {depth} using reaction checker")
                    nitro_reduction_found = True
                    nitro_reduction_depth = depth
                    reduced_nitro_mol = product
                else:
                    # Fallback to functional group checking
                    for reactant in reactants:
                        if checker.check_fg("Nitro group", reactant):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and checker.check_fg("Primary amine", product):
                                print(f"Nitro reduction detected at depth {depth} using FG checker")
                                nitro_reduction_found = True
                                nitro_reduction_depth = depth
                                reduced_nitro_mol = product

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal to find nitro reduction
    dfs_traverse(route)

    # If nitro reduction found, check for subsequent cyclization
    if nitro_reduction_found and nitro_reduction_depth is not None:
        print(f"Checking for cyclization after nitro reduction at depth {nitro_reduction_depth}")

        # Check reactions at lower depths (later in synthesis)
        for d, rxn_node in reaction_depths.items():
            if d < nitro_reduction_depth:  # Lower depth means later in synthesis
                rxn_rsmi = rxn_node["metadata"]["rsmi"]
                rxn_reactants = rxn_rsmi.split(">")[0].split(".")
                rxn_product = rxn_rsmi.split(">")[-1]

                # Check if this is a cyclization reaction
                cyclization_reactions = [
                    "Formation of NOS Heterocycles",
                    "Paal-Knorr pyrrole synthesis",
                    "Benzothiazole formation from aldehyde",
                    "Benzoxazole formation from aldehyde",
                    "Benzimidazole formation from aldehyde",
                    "Fischer indole",
                    "Friedlaender chinoline",
                ]

                for rxn_type in cyclization_reactions:
                    if checker.check_reaction(rxn_type, rxn_rsmi):
                        print(f"Cyclization reaction {rxn_type} detected at depth {d}")
                        subsequent_cyclization_found = True
                        break

                if not subsequent_cyclization_found:
                    # Fallback to ring counting
                    try:
                        reactants_rings = sum(
                            [
                                len(Chem.GetSSSR(Chem.MolFromSmiles(r)))
                                for r in rxn_reactants
                                if Chem.MolFromSmiles(r)
                            ]
                        )
                        product_mol = Chem.MolFromSmiles(rxn_product)
                        product_rings = len(Chem.GetSSSR(product_mol)) if product_mol else 0

                        if product_rings > reactants_rings:
                            print(
                                f"Cyclization detected at depth {d} by ring count: {reactants_rings} → {product_rings}"
                            )

                            # Check if any reactant contains the amine from nitro reduction
                            for reactant in rxn_reactants:
                                if reduced_nitro_mol and checker.check_fg(
                                    "Primary amine", reactant
                                ):
                                    print(f"Amine-containing reactant found in cyclization")
                                    subsequent_cyclization_found = True
                                    break
                    except Exception as e:
                        print(f"Error in ring counting: {e}")

    result = nitro_reduction_found and subsequent_cyclization_found
    print(
        f"Final result: {result} (Nitro reduction: {nitro_reduction_found}, Cyclization: {subsequent_cyclization_found})"
    )
    return result
