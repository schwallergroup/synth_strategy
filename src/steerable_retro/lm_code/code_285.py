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
    Detects if the synthesis involves a trifluoromethyl-substituted aromatic system
    """
    has_cf3_aromatic = False

    def dfs_traverse(node, depth=0):
        nonlocal has_cf3_aromatic

        if node["type"] == "mol" and node.get("smiles"):
            # Check for trifluoromethyl group attached to aromatic carbon using checker
            mol_smiles = node["smiles"]
            if checker.check_fg("Trifluoro group", mol_smiles) and any(
                checker.check_ring(ring, mol_smiles)
                for ring in [
                    "benzene",
                    "pyridine",
                    "naphthalene",
                    "anthracene",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                ]
            ):
                print(f"Depth {depth}: Detected trifluoromethyl aromatic in molecule: {mol_smiles}")
                has_cf3_aromatic = True

        elif node["type"] == "reaction" and node.get("metadata") and node["metadata"].get("rsmi"):
            # Check if the reaction introduces a trifluoromethyl group to an aromatic system
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has trifluoromethyl aromatic but reactants don't
            product_has_cf3 = checker.check_fg("Trifluoro group", product) and any(
                checker.check_ring(ring, product)
                for ring in [
                    "benzene",
                    "pyridine",
                    "naphthalene",
                    "anthracene",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                ]
            )

            if product_has_cf3:
                reactants_have_cf3 = any(
                    checker.check_fg("Trifluoro group", reactant)
                    and any(
                        checker.check_ring(ring, reactant)
                        for ring in [
                            "benzene",
                            "pyridine",
                            "naphthalene",
                            "anthracene",
                            "indole",
                            "quinoline",
                            "isoquinoline",
                        ]
                    )
                    for reactant in reactants
                    if reactant.strip()
                )

                if not reactants_have_cf3:
                    print(
                        f"Depth {depth}: Detected reaction introducing trifluoromethyl to aromatic: {rsmi}"
                    )
                    has_cf3_aromatic = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {has_cf3_aromatic}")

    return has_cf3_aromatic
