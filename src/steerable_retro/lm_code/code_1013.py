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
    Detects a convergent synthesis strategy where an indole derivative is coupled with a diarylmethanol
    in a late-stage reaction.
    """
    # Initialize tracking variables
    has_indole = False
    has_diarylmethanol = False
    has_late_coupling = False
    has_boc_protection = False

    def is_diarylmethanol(smiles):
        """Check if molecule is a diarylmethanol (secondary alcohol with two aromatic rings)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for secondary alcohol
        if not checker.check_fg("Secondary alcohol", smiles):
            return False

        # Find the carbon attached to OH and check if it's connected to two aromatic rings
        pattern = Chem.MolFromSmarts("[OH][C]")
        matches = mol.GetSubstructMatches(pattern)

        for match in matches:
            oh_idx, c_idx = match
            c_atom = mol.GetAtomWithIdx(c_idx)

            # Count aromatic ring connections to this carbon
            aromatic_ring_count = 0
            for neighbor in c_atom.GetNeighbors():
                if neighbor.GetIdx() != oh_idx and neighbor.GetIsAromatic():
                    aromatic_ring_count += 1

            if aromatic_ring_count >= 2:
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal has_indole, has_diarylmethanol, has_late_coupling, has_boc_protection

        if node["type"] == "mol":
            if "smiles" in node:
                # Check for indole substructure
                if checker.check_ring("indole", node["smiles"]):
                    has_indole = True
                    print(f"Found indole in molecule: {node['smiles']}")

                # Check for diarylmethanol
                if is_diarylmethanol(node["smiles"]):
                    has_diarylmethanol = True
                    print(f"Found diarylmethanol in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection
            if checker.check_reaction("Boc amine protection", rsmi) or checker.check_reaction(
                "Boc amine protection explicit", rsmi
            ):
                has_boc_protection = True
                print(f"Found Boc protection reaction: {rsmi}")

            # Check for late-stage coupling (depth 0 or 1)
            if depth <= 1:
                # Check if one reactant has indole and another has diarylmethanol
                has_indole_reactant = False
                has_diarylmethanol_reactant = False

                for reactant in reactants:
                    if checker.check_ring("indole", reactant):
                        has_indole_reactant = True
                    if is_diarylmethanol(reactant):
                        has_diarylmethanol_reactant = True

                # Check if product has both features
                product_has_indole = checker.check_ring("indole", product)
                product_has_diarylmethanol = is_diarylmethanol(product)

                # If reactants have separate features and product has both, it's a coupling
                if (
                    (has_indole_reactant and has_diarylmethanol_reactant)
                    or (has_indole_reactant and product_has_diarylmethanol)
                    or (has_diarylmethanol_reactant and product_has_indole)
                ):
                    has_late_coupling = True
                    print(f"Detected late-stage coupling of indole with diarylmethanol: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all required elements are present
    result = has_indole and has_diarylmethanol and has_late_coupling

    # Make Boc protection optional since it might not be necessary in all cases
    if has_boc_protection:
        print("Boc protection step found, which is typical for this strategy")
    else:
        print("No Boc protection found, but this might be an alternative approach")

    print(f"Indole detected: {has_indole}")
    print(f"Diarylmethanol detected: {has_diarylmethanol}")
    print(f"Late-stage coupling detected: {has_late_coupling}")
    print(f"Boc protection detected: {has_boc_protection}")
    print(f"Overall strategy detected: {result}")

    return result
