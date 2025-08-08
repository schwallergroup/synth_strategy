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
    This function detects if the synthesis route involves a piperazine scaffold.
    """
    piperazine_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_detected

        indent = "  " * depth
        print(f"{indent}Examining node at depth {depth}, type: {node['type']}")

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                print(f"{indent}Checking molecule: {mol_smiles}")

                # Use the checker function to detect piperazine
                if checker.check_ring("piperazine", mol_smiles):
                    print(f"{indent}✓ Detected piperazine scaffold in molecule: {mol_smiles}")
                    piperazine_detected = True
                else:
                    # Try to parse the molecule to see if there are any issues
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        # Check for piperazine-like patterns in the SMILES string
                        if "N1CC2CNCC" in mol_smiles or "N1CC2CC(CN" in mol_smiles:
                            print(
                                f"{indent}✓ Detected piperazine-like pattern in SMILES: {mol_smiles}"
                            )
                            piperazine_detected = True
                        else:
                            # Additional check with basic piperazine SMARTS pattern
                            piperazine_pattern = Chem.MolFromSmarts("N1CCNCC1")
                            if mol.HasSubstructMatch(piperazine_pattern):
                                print(
                                    f"{indent}✓ Detected piperazine scaffold using basic SMARTS pattern in: {mol_smiles}"
                                )
                                piperazine_detected = True
                            else:
                                # Try with the official pattern from checker
                                piperazine_smarts = checker.get_ring_smiles("piperazine")
                                if piperazine_smarts:
                                    piperazine_pattern = Chem.MolFromSmarts(piperazine_smarts)
                                    if piperazine_pattern and mol.HasSubstructMatch(
                                        piperazine_pattern
                                    ):
                                        print(
                                            f"{indent}✓ Detected piperazine scaffold using checker SMARTS pattern in: {mol_smiles}"
                                        )
                                        piperazine_detected = True
            except Exception as e:
                print(f"{indent}Error processing molecule node: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                # Check if any reactant or product contains piperazine
                rsmi = node["metadata"]["rsmi"]
                print(f"{indent}Checking reaction: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check product for piperazine
                if checker.check_ring("piperazine", product):
                    print(f"{indent}✓ Detected piperazine scaffold in reaction product: {product}")
                    piperazine_detected = True
                else:
                    # Check for piperazine-like patterns in product
                    if "N1CC2CNCC" in product or "N1CC2CC(CN" in product:
                        print(
                            f"{indent}✓ Detected piperazine-like pattern in product SMILES: {product}"
                        )
                        piperazine_detected = True

                # Check reactants for piperazine
                for reactant in reactants:
                    if checker.check_ring("piperazine", reactant):
                        print(
                            f"{indent}✓ Detected piperazine scaffold in reaction reactant: {reactant}"
                        )
                        piperazine_detected = True
                    else:
                        # Check for piperazine-like patterns in reactants
                        if "N1CC2CNCC" in reactant or "N1CC2CC(CN" in reactant:
                            print(
                                f"{indent}✓ Detected piperazine-like pattern in reactant SMILES: {reactant}"
                            )
                            piperazine_detected = True
            except Exception as e:
                print(f"{indent}Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting traversal to detect piperazine scaffold...")
    dfs_traverse(route)
    print(f"Piperazine scaffold detected: {piperazine_detected}")
    return piperazine_detected
