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
    Detects Friedel-Crafts acylation (aromatic C-H to C-C(=O)R) in the synthesis route.
    """
    acylation_found = False

    def dfs_traverse(node):
        nonlocal acylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # First, check if this is a Friedel-Crafts acylation reaction using the checker
            if checker.check_reaction("Friedel-Crafts acylation", rsmi):
                print(f"Friedel-Crafts acylation detected via reaction checker: {rsmi}")
                acylation_found = True
            else:
                # If the direct check fails, perform a more detailed analysis
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has a ketone group
                if checker.check_fg("Ketone", product):
                    print(f"Product contains ketone: {product}")

                    # Check if any reactant is aromatic
                    has_aromatic_reactant = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                for atom in mol.GetAtoms():
                                    if atom.GetIsAromatic():
                                        has_aromatic_reactant = True
                                        print(f"Found aromatic reactant: {reactant}")
                                        break
                        except:
                            continue

                    # Check if any reactant already has a ketone
                    ketone_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Ketone", reactant):
                            ketone_in_reactants = True
                            print(f"Ketone found in reactant: {reactant}")
                            break

                    # If we have an aromatic reactant, a ketone in the product, but no ketone in reactants,
                    # this is likely a Friedel-Crafts acylation
                    if has_aromatic_reactant and not ketone_in_reactants:
                        print(f"Friedel-Crafts acylation detected via pattern analysis: {rsmi}")
                        acylation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return acylation_found
