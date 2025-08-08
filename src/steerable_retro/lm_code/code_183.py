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
    This function detects the formation of a thiophene heterocycle from acyclic precursors.
    """
    heterocycle_formation_detected = False

    def dfs_traverse(node):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if thiophene or thiazole is in product
                product_has_thiophene = checker.check_ring("thiophene", product_smiles)
                product_has_thiazole = checker.check_ring("thiazole", product_smiles)
                product_has_benzothiophene = checker.check_ring("benzothiophene", product_smiles)

                if product_has_thiophene or product_has_thiazole or product_has_benzothiophene:
                    print(f"Product contains sulfur heterocycle: {product_smiles}")

                    # Check if reactants are acyclic (don't contain the same heterocycle)
                    reactants_have_thiophene = any(
                        checker.check_ring("thiophene", reactant) for reactant in reactants_smiles
                    )
                    reactants_have_thiazole = any(
                        checker.check_ring("thiazole", reactant) for reactant in reactants_smiles
                    )
                    reactants_have_benzothiophene = any(
                        checker.check_ring("benzothiophene", reactant)
                        for reactant in reactants_smiles
                    )

                    # If product has thiophene but reactants don't, or product has thiazole but reactants don't
                    if (
                        (product_has_thiophene and not reactants_have_thiophene)
                        or (product_has_thiazole and not reactants_have_thiazole)
                        or (product_has_benzothiophene and not reactants_have_benzothiophene)
                    ):
                        print("Reactants do not contain the same sulfur heterocycle")

                        # Check if this is a heterocycle formation reaction
                        # Look for common heterocycle formation reactions
                        is_heterocycle_formation = any(
                            [
                                checker.check_reaction("benzothiophene", rsmi),
                                checker.check_reaction("thiazole", rsmi),
                                checker.check_reaction("Paal-Knorr pyrrole", rsmi),
                                checker.check_reaction("Formation of NOS Heterocycles", rsmi),
                                checker.check_reaction("{benzothiazole}", rsmi),
                                checker.check_reaction("{benzothiophene}", rsmi),
                                checker.check_reaction("{thiazole}", rsmi),
                            ]
                        )

                        if is_heterocycle_formation:
                            print("Confirmed heterocycle formation reaction")
                            heterocycle_formation_detected = True
                        else:
                            # Check if reactants have sulfur-containing functional groups
                            # that could lead to thiophene/thiazole formation
                            has_sulfur_precursors = any(
                                checker.check_fg("Aliphatic thiol", reactant)
                                or checker.check_fg("Aromatic thiol", reactant)
                                or checker.check_fg("Monosulfide", reactant)
                                or checker.check_fg("Thiocarbonyl", reactant)
                                or checker.check_fg("Thiocyanate", reactant)
                                or checker.check_fg("Isothiocyanate", reactant)
                                or checker.check_fg("Thioamide", reactant)
                                or checker.check_fg("Thiourea", reactant)
                                or "S=" in reactant  # Thione groups
                                or "SC" in reactant  # General sulfur-carbon bonds
                                for reactant in reactants_smiles
                            )

                            if has_sulfur_precursors:
                                print(
                                    "Detected sulfur heterocycle formation from sulfur-containing precursors"
                                )
                                heterocycle_formation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_formation_detected
