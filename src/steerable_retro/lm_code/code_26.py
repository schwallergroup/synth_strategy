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
    This function detects a linear synthesis strategy maintaining a dichlorobenzene core.
    """
    # Track if we found the key features
    dichlorobenzene_core_maintained = True
    linear_synthesis = True
    step_count = 0

    def has_dichlorobenzene_core(smiles):
        """Helper function to check if a molecule has a dichlorobenzene core"""
        # Check for aromatic halide
        if not checker.check_fg("Aromatic halide", smiles):
            return False

        # Count chlorines in the molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Count chlorines
        cl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl")

        # Check if there's a benzene ring
        has_benzene = checker.check_ring("benzene", smiles)

        # Return true if we have a benzene ring and at least 2 chlorines
        return has_benzene and cl_count >= 2

    def dfs_traverse(node, depth=0):
        nonlocal dichlorobenzene_core_maintained, linear_synthesis, step_count

        if node["type"] == "reaction":
            step_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a linear step (one main reactant plus reagents)
            # Count significant reactants (those that don't start with '[' which are often reagents)
            significant_reactants = [r for r in reactants_smiles if not r.startswith("[")]
            if len(significant_reactants) > 1:
                linear_synthesis = False
                print(
                    f"Found non-linear step at depth {depth}: {len(significant_reactants)} significant reactants"
                )

            # Check for dichlorobenzene core in both reactants and products
            product_has_core = has_dichlorobenzene_core(product_smiles)
            reactant_has_core = any(has_dichlorobenzene_core(r) for r in reactants_smiles)

            if not (product_has_core and reactant_has_core):
                dichlorobenzene_core_maintained = False
                print(f"Dichlorobenzene core not maintained at depth {depth}")
                print(f"  Product has core: {product_has_core}")
                print(f"  Any reactant has core: {reactant_has_core}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if core is maintained throughout a linear synthesis with multiple steps
    return dichlorobenzene_core_maintained and linear_synthesis and step_count >= 3
