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
    This function detects if the synthesis route includes a late-stage nitration
    of an aromatic amine (conversion of -NH2 to -NO2 group).

    Note: In retrosynthesis, this appears as a reduction of nitro to amine.
    """
    nitration_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitration_found

        if node["type"] == "reaction" and depth <= 2:  # Late stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product_part = rsmi.split(">")[-1]
                products = product_part.split(".")

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # In retrosynthesis, we're looking for nitro reduction (nitro in reactants, amine in products)
                # This corresponds to nitration in forward synthesis

                # Check if any reactant has a nitro group
                nitro_reactant = None
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant):
                        print(f"Nitro group found in reactant: {reactant}")
                        nitro_reactant = reactant
                        break

                if nitro_reactant:
                    # Check if any product has an aromatic amine
                    for product in products:
                        if checker.check_fg("Aniline", product) or checker.check_fg(
                            "Primary amine", product
                        ):
                            print(f"Amine found in product: {product}")

                            # Verify this is a reduction reaction (nitro to amine)
                            # Check if the reaction is consistent with nitro reduction
                            # In retrosynthesis, this would be the reverse of nitration

                            # Get atom mapping of nitro group in reactant
                            nitro_mol = Chem.MolFromSmiles(nitro_reactant)
                            product_mol = Chem.MolFromSmiles(product)

                            if nitro_mol and product_mol:
                                # Check if this is a valid nitration/reduction by examining the reaction type
                                # or by checking if the reaction involves conversion between nitro and amine groups

                                # Since we're in retrosynthesis, we're looking for reduction reactions
                                # that would correspond to nitration in forward synthesis
                                if checker.check_reaction(
                                    "Reduction of nitro groups to amines", rsmi
                                ) or (
                                    checker.check_fg("Nitro group", nitro_reactant)
                                    and checker.check_fg("Aniline", product)
                                    or checker.check_fg("Primary amine", product)
                                ):

                                    print(
                                        f"Late-stage nitration of aromatic amine detected at depth {depth}"
                                    )
                                    nitration_found = True
                                    return  # Exit early once found

        # Continue DFS traversal
        for child in node.get("children", []):
            if not nitration_found:  # Only continue if we haven't found a match yet
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {nitration_found}")
    return nitration_found
