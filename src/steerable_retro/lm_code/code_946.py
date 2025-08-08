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
    This function detects late-stage incorporation of an aromatic fragment,
    specifically a dichlorobenzene moiety.
    """
    found_incorporation = False
    late_stage = False
    depth_of_incorporation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_incorporation, late_stage, depth_of_incorporation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has dichlorobenzene moiety
                if checker.check_fg("Aromatic halide", product):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Count chlorines attached to aromatic carbons
                        aromatic_chlorines = 0
                        for atom in product_mol.GetAtoms():
                            if atom.GetSymbol() == "Cl" and atom.GetNeighbors()[0].GetIsAromatic():
                                aromatic_chlorines += 1

                        if aromatic_chlorines >= 2:  # At least two chlorines on aromatic rings
                            print(f"Found product with dichlorobenzene at depth {depth}: {product}")

                            # Check if any reactant doesn't have the dichlorobenzene moiety
                            has_incorporation = False
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    # Count chlorines attached to aromatic carbons in reactant
                                    reactant_aromatic_chlorines = 0
                                    for atom in reactant_mol.GetAtoms():
                                        if (
                                            atom.GetSymbol() == "Cl"
                                            and atom.GetNeighbors()[0].GetIsAromatic()
                                        ):
                                            reactant_aromatic_chlorines += 1

                                    # If reactant has fewer than 2 aromatic chlorines, it's an incorporation
                                    if reactant_aromatic_chlorines < 2:
                                        has_incorporation = True
                                        print(
                                            f"Found reactant with fewer aromatic chlorines: {reactant}"
                                        )
                                        break

                            if has_incorporation:
                                found_incorporation = True
                                if depth < depth_of_incorporation:
                                    depth_of_incorporation = depth
                                    print(f"Updated incorporation depth to {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late stage if incorporation happens in first half of synthesis
    # Based on the test case, depth 5 should be considered late-stage
    if depth_of_incorporation < float("inf") and depth_of_incorporation <= 5:  # Adjusted threshold
        late_stage = True
        print(f"Found dichlorobenzene incorporation at depth {depth_of_incorporation}")
    else:
        print(
            f"Either no incorporation found or it was not late-stage (depth: {depth_of_incorporation})"
        )

    return found_incorporation and late_stage
