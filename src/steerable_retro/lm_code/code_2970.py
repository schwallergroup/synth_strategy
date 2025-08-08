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
    This function detects if aromatic groups are incorporated late in the synthesis
    (in the final step but not in earlier steps).
    """
    aromatic_incorporation_depths = []

    # Common aromatic rings to check
    aromatic_rings = [
        "benzene",
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "pyrazole",
        "oxazole",
        "thiazole",
        "naphthalene",
        "indole",
        "quinoline",
    ]

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains aromatic rings
                product_has_aromatic = False
                for ring in aromatic_rings:
                    if checker.check_ring(ring, product):
                        product_has_aromatic = True
                        break

                # If no aromatic rings found with checker, use RDKit's aromaticity detection
                if not product_has_aromatic:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and any(atom.GetIsAromatic() for atom in product_mol.GetAtoms()):
                        product_has_aromatic = True

                # Only proceed if product has aromatic rings
                if product_has_aromatic:
                    # Check if any reactant lacks aromatic rings
                    reactants_all_aromatic = True
                    for reactant in reactants:
                        reactant_has_aromatic = False
                        for ring in aromatic_rings:
                            if checker.check_ring(ring, reactant):
                                reactant_has_aromatic = True
                                break

                        # If no aromatic rings found with checker, use RDKit's aromaticity detection
                        if not reactant_has_aromatic:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and any(
                                atom.GetIsAromatic() for atom in reactant_mol.GetAtoms()
                            ):
                                reactant_has_aromatic = True

                        if not reactant_has_aromatic:
                            reactants_all_aromatic = False
                            break

                    # If not all reactants have aromatic rings but product does,
                    # then aromatic group was incorporated in this step
                    if not reactants_all_aromatic:
                        print(f"Found aromatic incorporation at depth {depth}")
                        print(f"Reaction: {rsmi}")
                        aromatic_incorporation_depths.append(depth)
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if aromatic incorporation happens only at the lowest depth (latest stage)
    if not aromatic_incorporation_depths:
        return False

    min_depth = min(aromatic_incorporation_depths)
    return all(depth == min_depth for depth in aromatic_incorporation_depths)
