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
    Detects synthesis routes based on tetrahydroisoquinoline scaffold
    """
    found_tetrahydroisoquinoline = False
    tetrahydroisoquinoline_reactions = 0

    def is_tetrahydroisoquinoline(smiles):
        """Helper function to check if a molecule contains tetrahydroisoquinoline scaffold"""
        if not smiles:
            return False

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # First check: Direct isoquinoline check
        if checker.check_ring("isoquinoline", smiles):
            # Ensure it's tetrahydro (not fully aromatic)
            # Look for non-aromatic nitrogen in a ring
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "N" and atom.IsInRing() and not atom.GetIsAromatic():
                    return True

        # Second check: Look for fused piperidine and benzene rings
        if checker.check_ring("piperidine", smiles) and checker.check_ring("benzene", smiles):
            # Get ring indices
            piperidine_indices = checker.get_ring_atom_indices("piperidine", smiles)
            benzene_indices = checker.get_ring_atom_indices("benzene", smiles)

            if piperidine_indices and benzene_indices:
                # Check if the rings share atoms (are fused)
                piperidine_atoms = set([atom for ring in piperidine_indices for atom in ring])
                benzene_atoms = set([atom for ring in benzene_indices for atom in ring])

                # If rings share atoms and there's a nitrogen in the piperidine ring
                if piperidine_atoms.intersection(benzene_atoms):
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "N" and atom.IsInRing() and not atom.GetIsAromatic():
                            return True

        # Third check: Direct pattern matching for 1,2,3,4-tetrahydroisoquinoline
        # This is a fallback in case the ring detection methods fail
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "N" and atom.IsInRing():
                # Check if it's part of a partially saturated isoquinoline-like structure
                if not atom.GetIsAromatic() and atom.GetDegree() in [2, 3]:
                    # This is a potential THIQ nitrogen
                    return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_tetrahydroisoquinoline, tetrahydroisoquinoline_reactions

        if node["type"] == "mol":
            # Check if the molecule contains tetrahydroisoquinoline scaffold
            if is_tetrahydroisoquinoline(node["smiles"]):
                found_tetrahydroisoquinoline = True
                print(f"Found tetrahydroisoquinoline scaffold at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                # Check if this reaction forms or modifies tetrahydroisoquinoline
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Handle multiple reactants
                reactants = reactants_part.split(".")

                # Check if product contains tetrahydroisoquinoline
                product_has_thiq = is_tetrahydroisoquinoline(product_part)

                # Check if any reactant contains tetrahydroisoquinoline
                reactants_have_thiq = any(
                    is_tetrahydroisoquinoline(reactant) for reactant in reactants
                )

                # Check for Pictet-Spengler reaction which commonly forms THIQ
                is_pictet_spengler = checker.check_reaction("Pictet-Spengler", rsmi)

                # If product has tetrahydroisoquinoline but reactants don't, this reaction forms it
                if product_has_thiq and not reactants_have_thiq:
                    tetrahydroisoquinoline_reactions += 1
                    reaction_type = "Pictet-Spengler" if is_pictet_spengler else "forming"
                    print(
                        f"Found reaction {reaction_type} tetrahydroisoquinoline at depth {depth}: {rsmi}"
                    )

                # If both product and reactants have tetrahydroisoquinoline, this reaction modifies it
                elif product_has_thiq and reactants_have_thiq:
                    tetrahydroisoquinoline_reactions += 1
                    print(
                        f"Found reaction modifying tetrahydroisoquinoline at depth {depth}: {rsmi}"
                    )
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found tetrahydroisoquinoline and at least one reaction involving it
    result = found_tetrahydroisoquinoline and tetrahydroisoquinoline_reactions > 0
    print(
        f"Final result: {result} (found_scaffold: {found_tetrahydroisoquinoline}, reactions: {tetrahydroisoquinoline_reactions})"
    )
    return result
