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
    Detects a strategy involving thiazole formation from a bromoacetyl intermediate
    and thioacetamide.
    """
    bromoacetyl_present = False
    thiazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal bromoacetyl_present, thiazole_formation

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for bromoacetyl group in reactants
                for reactant in reactants:
                    # Look for a primary halide that has a bromoacetyl pattern
                    if checker.check_fg("Primary halide", reactant):
                        mol = Chem.MolFromSmiles(reactant)
                        if (
                            mol
                            and mol.GetSubstructMatch(Chem.MolFromSmarts("[#6](=[O])[CH2]Br")) != ()
                        ):
                            bromoacetyl_present = True
                            print(f"Detected bromoacetyl group at depth {depth}")

                # Check for thioacetamide in reactants
                thioacetamide_present = any(
                    checker.check_fg("Thioamide", reactant) for reactant in reactants
                )

                # Check for thiazole in product
                thiazole_in_product = checker.check_ring("thiazole", product)

                # Check if this is a thiazole formation reaction
                thiazole_reaction = checker.check_reaction("thiazole", rsmi)

                if bromoacetyl_present and thioacetamide_present and thiazole_in_product:
                    thiazole_formation = True
                    print(
                        f"Detected thiazole formation from bromoacetyl and thioacetamide at depth {depth}"
                    )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if bromoacetyl_present and thiazole_formation:
        print("Detected bromoacetyl-thiazole formation strategy")
        return True

    return False
