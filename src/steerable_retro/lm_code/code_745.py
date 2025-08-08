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
    This function detects if the synthesis route involves early formation of a diaryl ether motif.
    """
    diaryl_ether_formed = False
    high_depth_threshold = 5  # Consider depths >= 5 as early in synthesis

    def is_diaryl_ether(mol_smiles):
        """Check if molecule contains a diaryl ether motif"""
        if not checker.check_fg("Ether", mol_smiles):
            return False

        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Find all ether bonds
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()

                # Check if one atom is oxygen and connected to two aromatic carbons
                if begin_atom.GetAtomicNum() == 8 and begin_atom.GetDegree() == 2:
                    neighbors = [mol.GetAtomWithIdx(n.GetIdx()) for n in begin_atom.GetNeighbors()]
                    if all(n.GetIsAromatic() and n.GetAtomicNum() == 6 for n in neighbors):
                        return True

                elif end_atom.GetAtomicNum() == 8 and end_atom.GetDegree() == 2:
                    neighbors = [mol.GetAtomWithIdx(n.GetIdx()) for n in end_atom.GetNeighbors()]
                    if all(n.GetIsAromatic() and n.GetAtomicNum() == 6 for n in neighbors):
                        return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal diaryl_ether_formed

        if node["type"] == "reaction":
            # Get depth from metadata or use traversal depth
            node_depth = int(node.get("metadata", {}).get("depth", depth))
            print(f"Checking reaction at depth {node_depth}")

            if node_depth >= high_depth_threshold:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reaction SMILES: {rsmi}")

                    # Check if this is a known ether-forming reaction
                    is_ether_reaction = (
                        checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction("Ullmann condensation", rsmi)
                        or checker.check_reaction(
                            "Ullmann-Goldberg Substitution aryl alcohol", rsmi
                        )
                    )

                    if is_ether_reaction:
                        print(f"Detected ether-forming reaction at depth {node_depth}")

                    # Check if product contains diaryl ether but reactants don't
                    product_has_diaryl_ether = is_diaryl_ether(product)

                    if product_has_diaryl_ether:
                        reactant_has_diaryl_ether = False
                        for reactant in reactants:
                            if is_diaryl_ether(reactant):
                                reactant_has_diaryl_ether = True
                                break

                        if not reactant_has_diaryl_ether:
                            print(f"Diaryl ether formation detected at depth {node_depth}")
                            diaryl_ether_formed = True

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return diaryl_ether_formed
