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
    Detects a convergent synthesis where multiple complex fragments
    are combined in the final or late-stage steps.
    """
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent

        # Check reaction nodes at late stages (final or near-final steps)
        if node["type"] == "reaction" and depth <= 2:
            # Extract reactants
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Skip if only one reactant (not convergent)
            if len(reactants) < 2:
                print(f"Depth {depth}: Not convergent - only one reactant")
                for child in node.get("children", []):
                    dfs_traverse(child, depth + 1)
                return

            # Check if we have complex reactants
            complex_reactants = 0
            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol:
                    # Define complexity based on atom count, rings, or functional groups
                    atom_count = mol.GetNumAtoms()
                    ring_count = mol.GetRingInfo().NumRings()

                    # Consider a reactant complex if it has >8 atoms or contains rings
                    if atom_count > 8 or ring_count > 0:
                        complex_reactants += 1
                        print(
                            f"Depth {depth}: Complex reactant found - atoms: {atom_count}, rings: {ring_count}"
                        )

            # Check for coupling reactions which are often used in convergent synthesis
            is_coupling = False
            rxn_smiles = node["metadata"].get("smiles", "")

            if rxn_smiles:
                coupling_reactions = [
                    "Suzuki",
                    "Negishi",
                    "Stille",
                    "Heck",
                    "Sonogashira",
                    "Buchwald-Hartwig",
                    "Ullmann",
                ]

                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rxn_smiles):
                        print(f"Depth {depth}: Coupling reaction detected: {rxn_type}")
                        is_coupling = True
                        break

            # Determine if this is a convergent step
            if complex_reactants >= 2:
                print(
                    f"Depth {depth}: Convergent synthesis detected with {complex_reactants} complex fragments"
                )
                is_convergent = True
            elif complex_reactants >= 1 and is_coupling:
                print(f"Depth {depth}: Convergent synthesis detected with coupling reaction")
                is_convergent = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if not is_convergent:
        print("No convergent synthesis pattern detected")

    return is_convergent
