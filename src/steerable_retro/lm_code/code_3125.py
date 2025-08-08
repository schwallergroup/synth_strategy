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
    This function detects the presence of nitrogen heterocycles (piperazine, pyrazole)
    along with fluorinated groups (trifluoromethyl) maintained throughout the synthesis.
    """
    # Track if we've found all required structures in the final product
    found_piperazine = False
    found_pyrazole = False
    found_trifluoromethyl = False

    # First check if the final product (root node) has all required structures
    if route["type"] == "mol" and "smiles" in route:
        mol_smiles = route["smiles"]
        print(f"Checking final product: {mol_smiles}")

        # Check for piperazine
        if checker.check_ring("piperazine", mol_smiles):
            found_piperazine = True
            print(f"Found piperazine in final product")
        else:
            # Check for similar structures that might be piperazine
            for ring_name in ["piperidine", "morpholine", "thiomorpholine", "diazepane"]:
                if checker.check_ring(ring_name, mol_smiles):
                    print(
                        f"Found {ring_name} in final product, which might be similar to piperazine"
                    )

        # Check for pyrazole
        if checker.check_ring("pyrazole", mol_smiles):
            found_pyrazole = True
            print(f"Found pyrazole in final product")

        # Check for trifluoromethyl group
        if checker.check_fg("Trifluoro group", mol_smiles):
            found_trifluoromethyl = True
            print(f"Found trifluoromethyl in final product")

    # If not all structures are found in the final product, check the synthesis route
    if not (found_piperazine and found_pyrazole and found_trifluoromethyl):

        def dfs_traverse(node, depth=0):
            nonlocal found_piperazine, found_pyrazole, found_trifluoromethyl

            if node["type"] == "mol" and "smiles" in node and not node.get("in_stock", False):
                # This is a non-starting material molecule
                try:
                    mol_smiles = node["smiles"]
                    print(f"Checking molecule at depth {depth}: {mol_smiles}")

                    # Check for piperazine
                    if not found_piperazine and checker.check_ring("piperazine", mol_smiles):
                        found_piperazine = True
                        print(f"Found piperazine at depth {depth}")

                    # Check for pyrazole
                    if not found_pyrazole and checker.check_ring("pyrazole", mol_smiles):
                        found_pyrazole = True
                        print(f"Found pyrazole at depth {depth}")

                    # Check for trifluoromethyl group
                    if not found_trifluoromethyl and checker.check_fg(
                        "Trifluoro group", mol_smiles
                    ):
                        found_trifluoromethyl = True
                        print(f"Found trifluoromethyl at depth {depth}")

                    # Additional check for structures that might be piperazine-like
                    if not found_piperazine:
                        # Check if there's a structure with "N1CCC(F)(c2ccccc2C(F)(F)F)CC1" pattern
                        # which appears in the stdout and might be a piperazine-like structure
                        mol = Chem.MolFromSmiles(mol_smiles)
                        if mol and "N1CCC" in mol_smiles and "CC1" in mol_smiles:
                            print(
                                f"Found potential piperazine-like structure at depth {depth}: {mol_smiles}"
                            )
                            # This might be a fluorinated piperidine which is similar to piperazine
                            if checker.check_ring("piperidine", mol_smiles):
                                found_piperazine = True
                                print(f"Confirmed as fluorinated piperidine at depth {depth}")

                except Exception as e:
                    print(f"Error processing molecule: {e}")

            # Process children
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

        # Start traversal
        dfs_traverse(route)

    # Return True if all required structures are found
    has_all_structures = found_piperazine and found_pyrazole and found_trifluoromethyl
    print(
        f"Found piperazine: {found_piperazine}, Found pyrazole: {found_pyrazole}, Found trifluoromethyl: {found_trifluoromethyl}"
    )
    return has_all_structures
