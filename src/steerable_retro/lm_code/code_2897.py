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
    This function detects if the synthesis route includes N-dealkylation
    (breaking of N-C bond in tertiary amine).
    """
    n_dealkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_dealkylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                # In retrosynthesis, the products are on the left side of the arrow
                # and the reactants are on the right side
                products = rsmi.split(">")[0].split(".")
                reactants = rsmi.split(">")[-1].split(".")

                # Check for specific N-dealkylation reactions
                if checker.check_reaction("Hydrogenolysis of tertiary amines", rsmi):
                    print(f"Found hydrogenolysis of tertiary amines reaction at depth {depth}")
                    n_dealkylation_found = True
                    return

                # Check for tertiary amine in reactant (target molecule in retrosynthesis)
                for reactant in reactants:
                    if checker.check_fg("Tertiary amine", reactant):
                        # Check if any product (starting material in retrosynthesis) has primary or secondary amine
                        # and the tertiary amine is no longer present in that product
                        for product in products:
                            if (
                                checker.check_fg("Primary amine", product)
                                or checker.check_fg("Secondary amine", product)
                            ) and not checker.check_fg("Tertiary amine", product):
                                print(f"N-dealkylation detected at depth {depth}")
                                print(f"Target molecule: {reactant}")
                                print(f"Starting material: {product}")
                                n_dealkylation_found = True
                                return

                # Check for other reactions that might involve N-dealkylation
                # For example, oxidative cleavage of tertiary amines
                for reactant in reactants:
                    if checker.check_fg("Tertiary amine", reactant):
                        for product in products:
                            # Check for secondary amine formation
                            if checker.check_fg(
                                "Secondary amine", product
                            ) and not checker.check_fg("Tertiary amine", product):
                                print(
                                    f"Potential N-dealkylation via secondary amine formation at depth {depth}"
                                )
                                n_dealkylation_found = True
                                return

                            # Check for primary amine formation
                            if checker.check_fg("Primary amine", product) and not checker.check_fg(
                                "Tertiary amine", product
                            ):
                                print(
                                    f"Potential N-dealkylation via primary amine formation at depth {depth}"
                                )
                                n_dealkylation_found = True
                                return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_dealkylation_found
