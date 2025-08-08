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
    This function detects if the synthesis preserves complex heterocyclic systems
    (trifluoromethyl pyrazole and fused thiophene-cyclohexane) throughout.
    """
    # Track if heterocycles are present in final product
    final_product_has_tf_pyrazole = False
    final_product_has_thiophene_cyclohexane = False

    # Track if heterocycles are broken during synthesis
    tf_pyrazole_broken = False
    thiophene_cyclohexane_broken = False

    def are_rings_fused(mol_smiles, ring1, ring2):
        """Helper function to check if two rings are fused in a molecule"""
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                return False

            # Get ring info
            ring_info = mol.GetRingInfo()

            # Get atom indices for each ring
            ring1_indices_list = checker.get_ring_atom_indices(ring1, mol_smiles)
            ring2_indices_list = checker.get_ring_atom_indices(ring2, mol_smiles)

            if not ring1_indices_list or not ring2_indices_list:
                return False

            # Check if any pair of rings share at least two atoms (fusion)
            for ring1_indices in ring1_indices_list:
                for ring2_indices in ring2_indices_list:
                    # Convert to sets for intersection
                    ring1_set = set(ring1_indices[0])
                    ring2_set = set(ring2_indices[0])

                    # If they share at least 2 atoms, they're fused
                    if len(ring1_set.intersection(ring2_set)) >= 2:
                        return True

            return False
        except Exception as e:
            print(f"Error checking fused rings: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_tf_pyrazole, final_product_has_thiophene_cyclohexane
        nonlocal tf_pyrazole_broken, thiophene_cyclohexane_broken

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for heterocycles in the molecule
            has_tf_pyrazole = checker.check_fg(
                "Trifluoro group", mol_smiles
            ) and checker.check_ring("pyrazole", mol_smiles)
            has_thiophene_cyclohexane = (
                checker.check_ring("thiophene", mol_smiles)
                and checker.check_ring("cyclohexane", mol_smiles)
                and are_rings_fused(mol_smiles, "thiophene", "cyclohexane")
            )

            # For final product (depth 0), check if heterocycles are present
            if depth == 0:
                final_product_has_tf_pyrazole = has_tf_pyrazole
                final_product_has_thiophene_cyclohexane = has_thiophene_cyclohexane
                print(
                    f"Final product has trifluoromethyl pyrazole: {final_product_has_tf_pyrazole}"
                )
                print(
                    f"Final product has fused thiophene-cyclohexane: {final_product_has_thiophene_cyclohexane}"
                )

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Analyze reaction to check if heterocycles are preserved
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if heterocycles are present in product (forward direction)
                product_has_tf_pyrazole = checker.check_fg(
                    "Trifluoro group", product
                ) and checker.check_ring("pyrazole", product)
                product_has_thiophene_cyclohexane = (
                    checker.check_ring("thiophene", product)
                    and checker.check_ring("cyclohexane", product)
                    and are_rings_fused(product, "thiophene", "cyclohexane")
                )

                # Check if heterocycles are present in reactants (forward direction)
                reactants_have_tf_pyrazole = any(
                    checker.check_fg("Trifluoro group", r) and checker.check_ring("pyrazole", r)
                    for r in reactants
                )
                reactants_have_thiophene_cyclohexane = any(
                    checker.check_ring("thiophene", r)
                    and checker.check_ring("cyclohexane", r)
                    and are_rings_fused(r, "thiophene", "cyclohexane")
                    for r in reactants
                )

                # Check if heterocycles are formed during reaction (not preserved in retrosynthesis)
                # In retrosynthesis, we're going from product to reactants
                if product_has_tf_pyrazole and not reactants_have_tf_pyrazole:
                    print(
                        f"Trifluoromethyl pyrazole formed in reaction at depth {depth} (broken in retrosynthesis)"
                    )
                    tf_pyrazole_broken = True

                if product_has_thiophene_cyclohexane and not reactants_have_thiophene_cyclohexane:
                    print(
                        f"Fused thiophene-cyclohexane formed in reaction at depth {depth} (broken in retrosynthesis)"
                    )
                    thiophene_cyclohexane_broken = True

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if heterocycles are preserved throughout synthesis
    tf_pyrazole_preserved = final_product_has_tf_pyrazole and not tf_pyrazole_broken
    thiophene_cyclohexane_preserved = (
        final_product_has_thiophene_cyclohexane and not thiophene_cyclohexane_broken
    )

    print(f"Trifluoromethyl pyrazole preserved: {tf_pyrazole_preserved}")
    print(f"Fused thiophene-cyclohexane preserved: {thiophene_cyclohexane_preserved}")

    # Return true only if both heterocycles are preserved or if they weren't in the final product
    return (tf_pyrazole_preserved or not final_product_has_tf_pyrazole) and (
        thiophene_cyclohexane_preserved or not final_product_has_thiophene_cyclohexane
    )
