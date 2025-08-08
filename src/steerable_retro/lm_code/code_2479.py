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
    Detects a strategy with multiple amide coupling reactions.
    """
    amide_coupling_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_count

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Depth {depth} - Analyzing reaction: {rsmi}")

                # Check for various amide coupling reactions using the checker function
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with primary amine to imide",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Aminolysis of esters",
                    "Schotten-Baumann to ester",
                ]

                # First check if any of the specific reaction types match
                is_amide_coupling = False
                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        amide_coupling_count += 1
                        print(
                            f"Depth {depth} - Amide coupling detected ({reaction_type}), total count: {amide_coupling_count}"
                        )
                        is_amide_coupling = True
                        break

                # If no specific reaction type matched, try a more general approach
                if not is_amide_coupling:
                    # Extract reactants and product
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Depth {depth} - Reactants: {reactants_smiles}")
                    print(f"Depth {depth} - Product: {product_smiles}")

                    # Check for reactant functional groups involved in amide formation
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants_smiles
                    )
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                    has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants_smiles)

                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants_smiles
                    )

                    # Count amides in reactants and product
                    reactant_amides = []
                    for r in reactants_smiles:
                        if checker.check_fg("Primary amide", r):
                            reactant_amides.append(("Primary amide", r))
                        if checker.check_fg("Secondary amide", r):
                            reactant_amides.append(("Secondary amide", r))
                        if checker.check_fg("Tertiary amide", r):
                            reactant_amides.append(("Tertiary amide", r))

                    product_amides = []
                    if checker.check_fg("Primary amide", product_smiles):
                        product_amides.append(("Primary amide", product_smiles))
                        print(f"Depth {depth} - Product has primary amide")
                    if checker.check_fg("Secondary amide", product_smiles):
                        product_amides.append(("Secondary amide", product_smiles))
                        print(f"Depth {depth} - Product has secondary amide")
                    if checker.check_fg("Tertiary amide", product_smiles):
                        product_amides.append(("Tertiary amide", product_smiles))
                        print(f"Depth {depth} - Product has tertiary amide")

                    reactant_amide_count = len(reactant_amides)
                    product_amide_count = len(product_amides)

                    print(
                        f"Depth {depth} - Reactants have: Carboxylic acid: {has_carboxylic_acid}, Acyl halide: {has_acyl_halide}, "
                        f"Ester: {has_ester}, Anhydride: {has_anhydride}, Amine: {has_amine}"
                    )
                    print(
                        f"Depth {depth} - Amide count - Reactants: {reactant_amide_count}, Product: {product_amide_count}"
                    )

                    # If product has at least as many amides as reactants and reactants have appropriate functional groups
                    if product_amide_count >= reactant_amide_count and (
                        (has_carboxylic_acid and has_amine)
                        or (has_acyl_halide and has_amine)
                        or (has_ester and has_amine)
                        or (has_anhydride and has_amine)
                    ):
                        # Additional check: ensure a new amide bond is formed
                        # This is a heuristic - if we have the right reactants and the product has an amide,
                        # it's likely an amide coupling
                        if product_amide_count > 0:
                            amide_coupling_count += 1
                            print(
                                f"Depth {depth} - Amide coupling detected (general approach), total count: {amide_coupling_count}"
                            )

            except Exception as e:
                print(f"Depth {depth} - Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Strategy criteria: Has at least 2 amide couplings
    result = amide_coupling_count >= 2
    print(
        f"Has multiple amide coupling strategy: {result} (found {amide_coupling_count} amide couplings)"
    )
    return result
