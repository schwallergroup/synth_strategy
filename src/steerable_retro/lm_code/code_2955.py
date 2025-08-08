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
    This function detects the preservation of a dichlorophenyl motif
    throughout the synthesis.
    """
    # Track dichlorophenyl-containing nodes and their depths
    dichlorophenyl_nodes = []

    def dfs_traverse(node, depth=0):
        # Check if the current node contains a dichlorophenyl motif
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Aromatic halide", mol_smiles):
                # Check specifically for dichlorophenyl pattern
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Count chlorine atoms attached to aromatic rings
                    aromatic_chlorines = 0
                    for atom in mol.GetAtoms():
                        if atom.GetIsAromatic() and atom.GetSymbol() == "C":
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetSymbol() == "Cl":
                                    aromatic_chlorines += 1

                    if aromatic_chlorines >= 2:
                        # Found a dichlorophenyl motif
                        dichlorophenyl_nodes.append((depth, mol_smiles))
                        print(f"Found dichlorophenyl motif at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # For reaction nodes, check if the dichlorophenyl motif is preserved
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if both reactants and product contain dichlorophenyl
                reactant_has_dichlorophenyl = False
                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant):
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            aromatic_chlorines = 0
                            for atom in mol.GetAtoms():
                                if atom.GetIsAromatic() and atom.GetSymbol() == "C":
                                    for neighbor in atom.GetNeighbors():
                                        if neighbor.GetSymbol() == "Cl":
                                            aromatic_chlorines += 1

                            if aromatic_chlorines >= 2:
                                reactant_has_dichlorophenyl = True
                                break

                product_has_dichlorophenyl = False
                if checker.check_fg("Aromatic halide", product):
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        aromatic_chlorines = 0
                        for atom in mol.GetAtoms():
                            if atom.GetIsAromatic() and atom.GetSymbol() == "C":
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.GetSymbol() == "Cl":
                                        aromatic_chlorines += 1

                        if aromatic_chlorines >= 2:
                            product_has_dichlorophenyl = True

                # If both have dichlorophenyl, the motif is preserved in this reaction
                if reactant_has_dichlorophenyl and product_has_dichlorophenyl:
                    print(f"Dichlorophenyl motif preserved in reaction at depth {depth}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if dichlorophenyl motif is preserved throughout synthesis
    # We need at least 2 nodes with dichlorophenyl at different depths
    if len(dichlorophenyl_nodes) >= 2:
        # Sort by depth to check if motif appears at different stages
        dichlorophenyl_nodes.sort(key=lambda x: x[0])
        depth_diff = dichlorophenyl_nodes[-1][0] - dichlorophenyl_nodes[0][0]

        # If nodes are at different depths, the motif is preserved through synthesis
        if depth_diff > 0:
            print(
                f"Dichlorophenyl motif preserved across {len(dichlorophenyl_nodes)} steps with depth difference of {depth_diff}"
            )
            return True

    print(
        f"Dichlorophenyl motif not sufficiently preserved: found in {len(dichlorophenyl_nodes)} nodes"
    )
    return False
