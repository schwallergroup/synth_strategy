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
    This function detects synthesis involving quinolone heterocyclic systems.
    """
    quinolone_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal quinolone_detected

        if node["type"] == "mol":
            smiles = node["smiles"]

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f"Could not parse molecule: {smiles}")
                    return

                # Check for quinolone structure - a quinoline with a carbonyl group
                if checker.check_ring("quinoline", smiles) or checker.check_ring(
                    "isoquinoline", smiles
                ):
                    # Check for carbonyl group in the right position (quinolone)
                    quinolone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)[nH]c2")
                    isoquinolone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c(=O)c2")

                    if mol.HasSubstructMatch(quinolone_pattern):
                        print(f"Quinolone system detected in molecule: {smiles}")
                        quinolone_detected = True
                    elif mol.HasSubstructMatch(isoquinolone_pattern):
                        print(f"Isoquinolone system detected in molecule: {smiles}")
                        quinolone_detected = True

                    # Additional check for N-substituted quinolones
                    n_subst_quinolone = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)nc2")
                    n_subst_isoquinolone = Chem.MolFromSmarts("c1ccc2c(c1)nc(=O)c2")

                    if mol.HasSubstructMatch(n_subst_quinolone):
                        print(f"N-substituted quinolone detected in molecule: {smiles}")
                        quinolone_detected = True
                    elif mol.HasSubstructMatch(n_subst_isoquinolone):
                        print(f"N-substituted isoquinolone detected in molecule: {smiles}")
                        quinolone_detected = True

            except Exception as e:
                print(f"Error analyzing molecule: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a known quinoline/quinolone forming reaction
                if checker.check_reaction("Friedlaender chinoline", rsmi):
                    print(f"Friedlaender quinoline synthesis detected: {rsmi}")
                    quinolone_detected = True

                # Check if product contains quinolone but reactants don't
                product_has_quinolone = False
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    quinolone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)[nH]c2")
                    isoquinolone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c(=O)c2")
                    n_subst_quinolone = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)nc2")
                    n_subst_isoquinolone = Chem.MolFromSmarts("c1ccc2c(c1)nc(=O)c2")

                    if (
                        product_mol.HasSubstructMatch(quinolone_pattern)
                        or product_mol.HasSubstructMatch(isoquinolone_pattern)
                        or product_mol.HasSubstructMatch(n_subst_quinolone)
                        or product_mol.HasSubstructMatch(n_subst_isoquinolone)
                    ):
                        product_has_quinolone = True

                reactants_have_quinolone = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        quinolone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)[nH]c2")
                        isoquinolone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c(=O)c2")
                        n_subst_quinolone = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)nc2")
                        n_subst_isoquinolone = Chem.MolFromSmarts("c1ccc2c(c1)nc(=O)c2")

                        if (
                            reactant_mol.HasSubstructMatch(quinolone_pattern)
                            or reactant_mol.HasSubstructMatch(isoquinolone_pattern)
                            or reactant_mol.HasSubstructMatch(n_subst_quinolone)
                            or reactant_mol.HasSubstructMatch(n_subst_isoquinolone)
                        ):
                            reactants_have_quinolone = True
                            break

                # If quinolone is formed in this reaction
                if product_has_quinolone and not reactants_have_quinolone:
                    print(f"Quinolone formation detected in reaction: {rsmi}")
                    quinolone_detected = True

                # Check for ring formation that might create quinoline structure
                if not reactants_have_quinolone and (
                    checker.check_ring("quinoline", product)
                    or checker.check_ring("isoquinoline", product)
                ):
                    print(f"Quinoline/isoquinoline ring formation detected: {rsmi}")
                    quinolone_detected = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return quinolone_detected
