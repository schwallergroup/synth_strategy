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
    Detects if the synthesis route involves a multi-fluorinated aromatic system.
    """
    multi_fluorinated_found = False

    def dfs_traverse(node, depth=0):
        nonlocal multi_fluorinated_found

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            try:
                # Check if molecule contains a trifluoro group (which has multiple F atoms)
                if checker.check_fg("Trifluoro group", mol_smiles):
                    print(f"Found trifluoro group in molecule: {mol_smiles}")
                    multi_fluorinated_found = True
                    return

                # Check if molecule contains aromatic halide with fluorine
                if checker.check_fg("Aromatic halide", mol_smiles):
                    # Create molecule to count fluorines on aromatic rings
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol is not None:
                        # Count fluorine atoms attached to aromatic carbons
                        aromatic_f_pattern = Chem.MolFromSmarts("[c]-[F]")
                        if aromatic_f_pattern:
                            matches = mol.GetSubstructMatches(aromatic_f_pattern)
                            f_aromatic_count = len(matches)

                            if f_aromatic_count >= 2:
                                print(
                                    f"Found multi-fluorinated aromatic system with {f_aromatic_count} F atoms in: {mol_smiles}"
                                )
                                multi_fluorinated_found = True
                                return
            except Exception as e:
                print(f"Error processing SMILES in fluorinated aromatic detection: {str(e)}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return multi_fluorinated_found
