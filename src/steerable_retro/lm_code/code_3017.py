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
    This function detects if an indole scaffold is preserved throughout the synthesis.
    It checks if the target molecule contains an indole scaffold and if this scaffold
    is preserved in the key intermediates of the synthesis.
    """
    # Track if indole is present in the target molecule and preserved in key steps
    target_has_indole = False
    indole_preserved = True

    # Track molecules that contain indole
    molecules_with_indole = []

    def dfs_traverse(node, depth=0):
        nonlocal target_has_indole, indole_preserved

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_indole = checker.check_ring("indole", mol_smiles)

            # If this is the target molecule (depth 0)
            if depth == 0:
                target_has_indole = has_indole
                print(f"Target molecule: {mol_smiles}")
                print(f"Target has indole: {target_has_indole}")

                # If target doesn't have indole, no need to check preservation
                if not target_has_indole:
                    indole_preserved = False
                    return

                molecules_with_indole.append(mol_smiles)

            # For intermediates, track if they have indole
            elif has_indole:
                print(f"Intermediate with indole found: {mol_smiles}")
                molecules_with_indole.append(mol_smiles)

        elif node["type"] == "reaction" and target_has_indole:
            # Check if this reaction preserves the indole scaffold
            # For reactions where the product has indole, at least one reactant should have indole
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                product_has_indole = checker.check_ring("indole", product)

                if product_has_indole:
                    # Check if any reactant has indole
                    reactant_has_indole = any(checker.check_ring("indole", r) for r in reactants)

                    # If product has indole but no reactant has it, indole is being formed
                    if not reactant_has_indole:
                        print(f"Indole is being formed in reaction: {rsmi}")
                    else:
                        print(f"Indole is preserved in reaction: {rsmi}")
                else:
                    # If a previous intermediate had indole but product doesn't, indole is not preserved
                    for reactant in reactants:
                        if checker.check_ring("indole", reactant):
                            print(f"Indole is lost in reaction: {rsmi}")
                            indole_preserved = False
                            break
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is successful if target has indole and indole is preserved throughout synthesis
    strategy_successful = target_has_indole and indole_preserved
    print(f"Indole scaffold preservation strategy detected: {strategy_successful}")

    return strategy_successful
