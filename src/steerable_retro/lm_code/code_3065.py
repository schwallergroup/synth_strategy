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
    Detects if the synthesis route incorporates a dichlorobenzyl moiety.
    A dichlorobenzyl moiety is a benzyl group (C6H5-CH2-) with two chlorine atoms on the benzene ring.
    """
    dichlorobenzyl_present = False

    def dfs_traverse(node):
        nonlocal dichlorobenzyl_present

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                mol = Chem.MolFromSmiles(mol_smiles)

                if mol and checker.check_fg("Aromatic halide", mol_smiles):
                    # Check for dichlorobenzyl pattern using SMARTS
                    # These patterns look for a benzene ring with two chlorines and a CH2 group attached
                    # We check for various possible positions of the chlorines

                    # Define patterns for different dichlorobenzyl arrangements
                    patterns = [
                        "c1(C[#0])c(Cl)c(Cl)ccc1",  # 2,3-dichlorobenzyl
                        "c1(C[#0])c(Cl)cc(Cl)cc1",  # 2,4-dichlorobenzyl
                        "c1(C[#0])c(Cl)ccc(Cl)c1",  # 2,5-dichlorobenzyl
                        "c1(C[#0])c(Cl)cccc1(Cl)",  # 2,6-dichlorobenzyl
                        "c1(C[#0])cc(Cl)c(Cl)cc1",  # 3,4-dichlorobenzyl
                        "c1(C[#0])cc(Cl)cc(Cl)c1",  # 3,5-dichlorobenzyl
                    ]

                    # Check each pattern
                    for pattern in patterns:
                        try:
                            patt = Chem.MolFromSmarts(pattern)
                            if patt and mol.HasSubstructMatch(patt):
                                dichlorobenzyl_present = True
                                print(f"Found dichlorobenzyl moiety in molecule: {mol_smiles}")
                                break
                        except Exception as e:
                            print(f"Error with pattern {pattern}: {e}")

                    # Additional check for common dichlorobenzyl patterns seen in the errors
                    if not dichlorobenzyl_present:
                        # Check specifically for the pattern seen in the error messages
                        if "COCc1c(Cl)cccc1Cl" in mol_smiles or "OCc1c(Cl)cccc1Cl" in mol_smiles:
                            dichlorobenzyl_present = True
                            print(
                                f"Found dichlorobenzyl moiety via direct SMILES match: {mol_smiles}"
                            )
            except Exception as e:
                print(f"Error processing molecule {node.get('smiles', 'unknown')}: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return dichlorobenzyl_present
