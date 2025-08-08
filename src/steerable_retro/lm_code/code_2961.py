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
    This function detects a linear synthesis strategy where fragments are added sequentially
    without convergent steps (each reaction has only one non-reagent reactant).
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node.get("type") == "reaction":
            if "metadata" in node and "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Count substantial reactants (not small reagents)
                substantial_reactants = []
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if (
                        r_mol and r_mol.GetNumHeavyAtoms() > 7
                    ):  # Increased threshold for "substantial"
                        substantial_reactants.append(r)
                        print(
                            f"  Substantial reactant found: {r} with {r_mol.GetNumHeavyAtoms()} heavy atoms"
                        )

                # Check if this is a fragment assembly reaction
                if len(substantial_reactants) > 1:
                    # Check if this is a coupling reaction
                    if (
                        checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation", rsmi
                        )
                        or checker.check_reaction("Sonogashira", rsmi)
                        or checker.check_reaction("Stille reaction", rsmi)
                        or checker.check_reaction("Negishi coupling", rsmi)
                        or checker.check_reaction("Heck", rsmi)
                    ):
                        print(
                            f"Found convergent coupling step at depth {depth} with {len(substantial_reactants)} substantial reactants"
                        )
                        is_linear = False
                    # Check for other fragment assembly reactions
                    elif (
                        checker.check_reaction("Amide formation", rsmi)
                        or checker.check_reaction("Esterification", rsmi)
                        or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction("Urea synthesis", rsmi)
                        or checker.check_reaction("Reductive amination", rsmi)
                    ):
                        print(
                            f"Found convergent fragment assembly step at depth {depth} with {len(substantial_reactants)} substantial reactants"
                        )
                        is_linear = False
                    else:
                        # Generic check for any reaction with multiple substantial reactants
                        print(
                            f"Found potential convergent step at depth {depth} with {len(substantial_reactants)} substantial reactants"
                        )
                        is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear fragment assembly strategy: {is_linear}")
    return is_linear
