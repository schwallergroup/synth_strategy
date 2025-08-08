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
    This function detects if the synthesis is based on forming a pyrazolone core
    and preserves specific substituents (like dichlorophenyl) throughout.
    """
    # Track if we've found a pyrazolone formation reaction
    found_pyrazolone_formation = False
    # Track if dichlorophenyl is preserved throughout relevant reactions
    dichlorophenyl_preserved = False
    # Check if final product has both pyrazolone and dichlorophenyl
    final_has_both = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pyrazolone_formation, dichlorophenyl_preserved, final_has_both

        # Check final product (depth 0)
        if depth == 0 and node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_pyrazolone = checker.check_ring("pyrazole", mol_smiles)
            has_dichlorophenyl = checker.check_fg("Aromatic halide", mol_smiles)

            if has_pyrazolone and has_dichlorophenyl:
                print(f"Final product has both pyrazolone and dichlorophenyl: {mol_smiles}")
                final_has_both = True
            else:
                print(
                    f"Final product missing required structures. Pyrazolone: {has_pyrazolone}, Dichlorophenyl: {has_dichlorophenyl}"
                )
                return

        # Check reaction nodes for pyrazolone formation
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]
            reactants = rxn_smiles.split(">")[0]
            product = rxn_smiles.split(">")[-1]

            # Check if this reaction forms a pyrazolone
            product_has_pyrazolone = checker.check_ring("pyrazole", product)
            reactants_have_pyrazolone = any(
                checker.check_ring("pyrazole", r) for r in reactants.split(".")
            )

            if product_has_pyrazolone and not reactants_have_pyrazolone:
                print(f"Found pyrazolone formation reaction: {rxn_smiles}")
                found_pyrazolone_formation = True

            # Check if dichlorophenyl is preserved in this reaction
            if checker.check_fg("Aromatic halide", product):
                for reactant in reactants.split("."):
                    if checker.check_fg("Aromatic halide", reactant):
                        dichlorophenyl_preserved = True
                        print(f"Dichlorophenyl preserved in reaction: {rxn_smiles}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return true if all conditions are met
    result = final_has_both and found_pyrazolone_formation and dichlorophenyl_preserved
    print(
        f"Final result: {result} (final_has_both: {final_has_both}, found_pyrazolone_formation: {found_pyrazolone_formation}, dichlorophenyl_preserved: {dichlorophenyl_preserved})"
    )
    return result
