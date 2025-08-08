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
    This function detects if a synthetic route employs a late-stage nitro reduction strategy,
    where an aromatic nitro group is reduced to an amine in the final or penultimate step.
    """
    # Initialize tracking variables
    has_late_stage_nitro_reduction = False
    nitro_present_in_intermediates = False

    # List of aromatic rings to check
    aromatic_rings = [
        "benzene",
        "naphthalene",
        "anthracene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "quinoline",
        "isoquinoline",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_nitro_reduction, nitro_present_in_intermediates

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if this is a nitro reduction reaction
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                # Alternative check: reactants have nitro groups, product has amines, and nitro groups are gone
                if not is_nitro_reduction:
                    reactants_have_nitro = any(
                        checker.check_fg("Nitro group", r) for r in reactants_smiles
                    )
                    product_has_amine = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )
                    product_lost_nitro = not checker.check_fg("Nitro group", product_smiles)

                    is_nitro_reduction = (
                        reactants_have_nitro and product_has_amine and product_lost_nitro
                    )

                if is_nitro_reduction:
                    print(f"Found nitro reduction reaction at depth {depth}")

                    # Check if the nitro group is on an aromatic ring in reactants
                    aromatic_nitro = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Nitro group", reactant):
                            if any(checker.check_ring(ring, reactant) for ring in aromatic_rings):
                                aromatic_nitro = True
                                break

                    # Consider depth 0 or 1 as late-stage (final or penultimate step)
                    if depth <= 1 and aromatic_nitro:
                        print(f"This is a late-stage aromatic nitro reduction (depth {depth})")
                        has_late_stage_nitro_reduction = True
                else:
                    # Check for aromatic nitro groups in the product (intermediate compound)
                    if depth > 1 and checker.check_fg("Nitro group", product_smiles):
                        # Verify the nitro group is on an aromatic ring
                        if any(checker.check_ring(ring, product_smiles) for ring in aromatic_rings):
                            print(f"Found aromatic nitro group in intermediate at depth {depth}")
                            nitro_present_in_intermediates = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if we have a late-stage nitro reduction
    # and nitro groups were present in earlier intermediates
    strategy_present = has_late_stage_nitro_reduction and nitro_present_in_intermediates
    print(f"Late-stage nitro reduction strategy detected: {strategy_present}")
    return strategy_present
