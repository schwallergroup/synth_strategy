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
    This function detects if the synthesis involves no ring formations or openings.
    """
    # Track if we've found any ring changes
    no_ring_changes = True

    def dfs_traverse(node):
        nonlocal no_ring_changes

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            try:
                # Process reactants
                reactants_ring_count = 0
                reactants_mols = []
                for reactant in reactants_smiles.split("."):
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        reactants_mols.append(mol)
                        reactants_ring_count += mol.GetRingInfo().NumRings()

                # Process products
                products_ring_count = 0
                products_mols = []
                for product in products_smiles.split("."):
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        products_mols.append(mol)
                        products_ring_count += mol.GetRingInfo().NumRings()

                # Check if ring count changed
                if reactants_ring_count != products_ring_count:
                    no_ring_changes = False
                    print(
                        f"Found ring formation/opening: {reactants_ring_count} -> {products_ring_count} in reaction: {rsmi}"
                    )
                else:
                    # Even if ring count is the same, check if ring types changed
                    # This could indicate ring opening and formation in the same reaction
                    ring_types = [
                        "furan",
                        "pyran",
                        "benzene",
                        "pyridine",
                        "piperidine",
                        "cyclohexane",
                        "cyclopentane",
                        "indole",
                        "pyrrole",
                    ]

                    reactants_rings = set()
                    products_rings = set()

                    for ring_type in ring_types:
                        for reactant in reactants_smiles.split("."):
                            if checker.check_ring(ring_type, reactant):
                                reactants_rings.add(ring_type)

                        for product in products_smiles.split("."):
                            if checker.check_ring(ring_type, product):
                                products_rings.add(ring_type)

                    if reactants_rings != products_rings:
                        no_ring_changes = False
                        print(
                            f"Found ring type change: {reactants_rings} -> {products_rings} in reaction: {rsmi}"
                        )

            except Exception as e:
                print(f"Error processing reaction SMILES: {rsmi}, Error: {str(e)}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"No ring formations/openings: {no_ring_changes}")
    return no_ring_changes
