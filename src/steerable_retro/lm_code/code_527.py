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
    This function detects if the synthetic route involves a late-stage deprotection
    of a tert-butyl ester to form a carboxylic acid.
    """
    late_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal late_deprotection

        if node["type"] == "reaction" and depth <= 2:  # Check final and near-final reactions
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is a deprotection reaction
                is_deprotection = False

                # Try specific reaction types first
                if checker.check_reaction("COOH ethyl deprotection", rsmi):
                    print("Detected COOH ethyl deprotection")
                    is_deprotection = True
                elif checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                    print("Detected Ester saponification (alkyl deprotection)")
                    is_deprotection = True
                elif checker.check_reaction("Tert-butyl deprotection of amine", rsmi):
                    print("Detected Tert-butyl deprotection of amine")
                    is_deprotection = True
                elif checker.check_reaction("Carboxyl benzyl deprotection", rsmi):
                    print("Detected Carboxyl benzyl deprotection")
                    is_deprotection = True

                # If no specific reaction type matched, check for the transformation pattern
                if not is_deprotection:
                    # Check if reactant has ester and product has carboxylic acid
                    has_ester_reactant = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                    has_acid_product = checker.check_fg("Carboxylic acid", product_smiles)

                    if has_ester_reactant and has_acid_product:
                        print("Detected ester to carboxylic acid transformation")
                        is_deprotection = True

                if is_deprotection:
                    print(f"Deprotection reaction detected at depth {depth}")

                    # Check for tert-butyl ester in reactants
                    has_tbutyl_ester = False

                    # Define tert-butyl ester pattern
                    tbutyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC(C)(C)C")

                    for r in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(tbutyl_ester_pattern):
                            print(f"Found tert-butyl ester in reactant: {r}")
                            has_tbutyl_ester = True
                            break

                    # Check for carboxylic acid in product
                    has_carboxylic_acid = checker.check_fg("Carboxylic acid", product_smiles)

                    print(
                        f"Has tert-butyl ester: {has_tbutyl_ester}, Has carboxylic acid: {has_carboxylic_acid}"
                    )

                    if has_tbutyl_ester and has_carboxylic_acid and depth <= 2:
                        late_deprotection = True
                        print("Late-stage tert-butyl ester deprotection detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_deprotection
