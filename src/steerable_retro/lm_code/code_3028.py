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
    This function detects a synthetic strategy involving a sequence of nitrogen
    functional group interconversions (nitro → amine → secondary amine).
    """
    # Track functional group transformations
    transformations = []

    def dfs_traverse(node, path=None):
        if path is None:
            path = []

        current_path = path.copy()

        if node["type"] == "mol":
            current_path.append(node["smiles"])

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            depth = node.get("metadata", {}).get("depth", -1)

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # In retrosynthesis, product is the starting material and reactants are the results
            # Check for nitro reduction (nitro → amine)
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Detected nitro reduction reaction at depth {depth}")
                transformations.append(("nitro_to_amine", depth))

            # Check for primary amine alkylation (primary amine → secondary amine)
            # Include reductive amination which also converts primary to secondary amines
            elif checker.check_reaction(
                "N-alkylation of primary amines with alkyl halides", rsmi
            ) or checker.check_reaction("Reductive amination with aldehyde", rsmi):
                print(f"Detected primary amine alkylation/reductive amination at depth {depth}")
                transformations.append(("amine_to_secondary_amine", depth))

            # If specific reaction checks fail, fall back to functional group analysis
            else:
                # Create molecules
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and any(reactant_mols):
                    # Check functional groups using checker functions
                    product_has_nitro = checker.check_fg("Nitro group", product)
                    product_has_primary_amine = checker.check_fg("Primary amine", product)
                    product_has_secondary_amine = checker.check_fg("Secondary amine", product)

                    reactant_has_nitro = any(
                        checker.check_fg("Nitro group", r) for r in reactants if r
                    )
                    reactant_has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants if r
                    )
                    reactant_has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants if r
                    )

                    # In retrosynthesis: product (starting material) → reactants (results)
                    if product_has_primary_amine and reactant_has_nitro:
                        print(f"Detected amine to nitro transformation at depth {depth}")
                        transformations.append(("nitro_to_amine", depth))
                    elif product_has_secondary_amine and reactant_has_primary_amine:
                        print(
                            f"Detected secondary amine to primary amine transformation at depth {depth}"
                        )
                        transformations.append(("amine_to_secondary_amine", depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_path)

    # Start traversal from root
    dfs_traverse(route)

    # Sort transformations by depth (highest to lowest, as higher depth is earlier in synthesis)
    # If depth is -1, treat it as a very high value to place it at the beginning (earliest)
    transformations.sort(key=lambda x: -float("inf") if x[1] == -1 else -x[1])

    print(f"All transformations (sorted by depth): {transformations}")

    # Check if we have the nitro → amine → secondary amine sequence
    # In retrosynthesis, this would appear as secondary_amine → amine → nitro
    transformation_types = [t[0] for t in transformations]

    if (
        "nitro_to_amine" in transformation_types
        and "amine_to_secondary_amine" in transformation_types
    ):
        nitro_idx = transformation_types.index("nitro_to_amine")
        amine_idx = transformation_types.index("amine_to_secondary_amine")

        # Check if the transformations are in the correct order
        # In retrosynthesis, amine_to_secondary_amine should come before nitro_to_amine
        # This means amine_idx should be less than nitro_idx in the sorted list
        if amine_idx < nitro_idx:
            print(
                "Detected nitrogen functional group interconversion sequence: nitro → amine → secondary amine"
            )
            return True

    return False
