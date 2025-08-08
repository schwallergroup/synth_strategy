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
    Detects if the synthetic route involves thiazole ring formation from a chloroacetyl compound and thiourea.
    """
    thiazole_formed = False
    thiourea_used = False
    chloroacetyl_used = False

    def dfs_traverse(node):
        nonlocal thiazole_formed, thiourea_used, chloroacetyl_used

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thiourea in reactants
                for reactant in reactants:
                    try:
                        if checker.check_fg("Thiourea", reactant):
                            thiourea_used = True
                            print(f"Found thiourea in reactant: {reactant}")
                    except Exception as e:
                        print(f"Error checking thiourea: {e}")
                        continue

                # Check for chloroacetyl compound in reactants (has both primary halide and ketone)
                for reactant in reactants:
                    try:
                        # Look for a molecule with both a primary halide (specifically chloride) and a ketone/acyl group
                        mol = Chem.MolFromSmiles(reactant)
                        if (
                            mol
                            and "Cl" in reactant
                            and checker.check_fg("Primary halide", reactant)
                            and checker.check_fg("Ketone", reactant)
                        ):
                            chloroacetyl_used = True
                            print(f"Found chloroacetyl compound in reactant: {reactant}")
                    except Exception as e:
                        print(f"Error checking chloroacetyl: {e}")
                        continue

                # Check for thiazole formation (thiazole in product but not in reactants)
                try:
                    # First check if thiazole is in any reactants
                    thiazole_in_reactants = any(
                        checker.check_ring("thiazole", r) for r in reactants
                    )

                    # Then check if thiazole is in the product
                    if checker.check_ring("thiazole", product) and not thiazole_in_reactants:
                        thiazole_formed = True
                        print(f"Found thiazole formation in product: {product}")

                    # Alternative: directly check if this is a thiazole formation reaction
                    if checker.check_reaction("thiazole", rsmi):
                        thiazole_formed = True
                        print(f"Detected thiazole formation reaction: {rsmi}")
                except Exception as e:
                    print(f"Error checking thiazole formation: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if both patterns were found and thiazole was formed
    result = thiazole_formed and thiourea_used and chloroacetyl_used
    print(f"Thiazole formation from chloroacetyl detected: {result}")
    return result
