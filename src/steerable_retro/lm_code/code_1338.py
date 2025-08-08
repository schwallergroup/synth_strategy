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
    Detects synthesis routes that produce fluorine-rich final products.
    """
    # Initialize flag
    has_fluorine_rich_product = False

    # Process the root node (final product)
    if route["type"] == "mol" and route.get("smiles"):
        try:
            smiles = route["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Count fluorine atoms directly
                fluorine_count = 0
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() == 9:  # Atomic number of fluorine
                        fluorine_count += 1

                # Check for specific fluorine-containing groups
                has_cf3 = checker.check_fg("Trifluoro group", smiles)
                has_ccl3 = checker.check_fg("Trichloro group", smiles)

                # Check for CF2H groups (difluoromethyl)
                cf2h_pattern = Chem.MolFromSmarts("[#6]([F])([F])[H]")
                cf2h_count = len(mol.GetSubstructMatches(cf2h_pattern)) if cf2h_pattern else 0

                # Check for other fluorinated groups
                has_fluoroalkyl = False
                if fluorine_count > 0:
                    # Check if any carbon is bonded to fluorine
                    for atom in mol.GetAtoms():
                        if atom.GetAtomicNum() == 6:  # Carbon
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetAtomicNum() == 9:  # Fluorine
                                    has_fluoroalkyl = True
                                    break

                # Define criteria for "fluorine-rich"
                if (
                    fluorine_count >= 3
                    or has_cf3
                    or (fluorine_count > 0 and cf2h_count > 0)
                    or (has_fluoroalkyl and fluorine_count >= 2)
                ):
                    print(f"Found fluorine-rich product with {fluorine_count} F atoms")
                    if has_cf3:
                        print("Contains CF3 group")
                    if cf2h_count > 0:
                        print(f"Contains {cf2h_count} CF2H groups")
                    has_fluorine_rich_product = True

        except Exception as e:
            print(f"Error processing molecule SMILES for fluorine content detection: {e}")

    return has_fluorine_rich_product
