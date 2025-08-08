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
    This function detects a strategy where a ketone is temporarily converted
    to other functional groups and then returned to a ketone
    """
    # Track the sequence of functional group transformations
    # In synthesis direction: [(source_fg, target_fg, depth, atom_indices), ...]
    transformations = []

    # List of functional groups that ketones can be converted to
    intermediate_fgs = [
        "Substituted imine",
        "Unsubstituted imine",
        "Hydrazone",
        "Oxime",
        "Primary alcohol",
        "Secondary alcohol",
        "Alkene",
        "Acetal/Ketal",
        "Enol",
        "Enamine",
        "Thiocarbonyl",
    ]

    def get_ketone_atom_indices(mol_smiles):
        """Get atom indices of ketone carbonyls in a molecule"""
        if not mol_smiles:
            return []
        return checker.get_fg_atom_indices("Ketone", mol_smiles)

    def get_fg_atom_indices(fg_name, mol_smiles):
        """Get atom indices of a functional group in a molecule"""
        if not mol_smiles:
            return []
        return checker.get_fg_atom_indices(fg_name, mol_smiles)

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check product for ketone (retrosynthetic: target molecule)
            product_ketone_indices = get_ketone_atom_indices(product)
            product_has_ketone = len(product_ketone_indices) > 0

            # Check reactants for ketone (retrosynthetic: starting materials)
            reactant_ketone_indices = []
            for r in reactants:
                if r:
                    indices = get_ketone_atom_indices(r)
                    if indices:
                        reactant_ketone_indices.extend(indices)
            reactant_has_ketone = len(reactant_ketone_indices) > 0

            # Check for intermediate functional groups
            product_intermediate_fgs = []
            for fg in intermediate_fgs:
                if checker.check_fg(fg, product):
                    product_intermediate_fgs.append((fg, get_fg_atom_indices(fg, product)))

            reactant_intermediate_fgs = []
            for r in reactants:
                if r:
                    for fg in intermediate_fgs:
                        if checker.check_fg(fg, r):
                            reactant_intermediate_fgs.append((fg, get_fg_atom_indices(fg, r)))

            # Record transformations in synthesis direction
            if product_has_ketone and reactant_intermediate_fgs:
                # Ketone in product, intermediate FG in reactants (retrosynthetic)
                # In synthesis direction: intermediate FG -> ketone
                for fg, indices in reactant_intermediate_fgs:
                    transformations.append((fg, "ketone", depth, product_ketone_indices))
                    print(f"Detected transformation: {fg} to ketone at depth {depth}")

            elif product_intermediate_fgs and reactant_has_ketone:
                # Intermediate FG in product, ketone in reactants (retrosynthetic)
                # In synthesis direction: ketone -> intermediate FG
                for fg, indices in product_intermediate_fgs:
                    transformations.append(("ketone", fg, depth, reactant_ketone_indices))
                    print(f"Detected transformation: ketone to {fg} at depth {depth}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth to get chronological order in synthesis direction
    transformations.sort(key=lambda x: x[2], reverse=True)

    print(f"All transformations (sorted by depth): {transformations}")

    # Check if we have a ketone → other functional groups → ketone pattern
    if len(transformations) < 2:
        print("Not enough transformations found")
        return False

    # Find the first transformation from ketone to something else
    ketone_to_intermediate_idx = -1
    for i, (source, target, _, _) in enumerate(transformations):
        if source == "ketone" and target in intermediate_fgs:
            ketone_to_intermediate_idx = i
            break

    if ketone_to_intermediate_idx == -1:
        print("No transformation from ketone to intermediate found")
        return False

    # Find a subsequent transformation from something back to ketone
    intermediate_to_ketone_found = False
    for i in range(ketone_to_intermediate_idx + 1, len(transformations)):
        source, target, _, _ = transformations[i]
        if source in intermediate_fgs and target == "ketone":
            intermediate_to_ketone_found = True
            break

    result = intermediate_to_ketone_found
    print(f"Ketone to intermediate found: {ketone_to_intermediate_idx != -1}")
    print(f"Intermediate to ketone found: {intermediate_to_ketone_found}")
    print(f"Final result: {result}")

    return result
