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
    Detects if a tricyclic scaffold is preserved throughout the synthesis.
    """
    # First identify if the final product has a tricyclic scaffold
    final_product = None

    def find_final_product(node):
        nonlocal final_product
        if node["type"] == "mol" and final_product is None:
            final_product = node["smiles"]
            return
        for child in node.get("children", []):
            find_final_product(child)

    find_final_product(route)

    if not final_product:
        print("No final product found")
        return False

    final_mol = Chem.MolFromSmiles(final_product)
    if not final_mol:
        print(f"Could not parse final product SMILES: {final_product}")
        return False

    # Check if final product has a tricyclic scaffold
    ring_info = final_mol.GetRingInfo()
    print(f"Final product has {ring_info.NumRings()} rings")

    # Expanded list of tricyclic scaffolds to check
    tricyclic_rings = [
        "anthracene",
        "phenanthrene",
        "acridine",
        "carbazole",
        "dibenzofuran",
        "dibenzothiophene",
        "phenothiazine",
        "phenoxazine",
        "xanthene",
        "thioxanthene",
        "purine",
        "pteridin",
        "benzotriazole",
        "indazole",
        "benzimidazole",
        "benzothiazole",
        "benzoxazole",
    ]

    # Check if final product contains any of these tricyclic scaffolds
    final_product_scaffold = None
    for ring in tricyclic_rings:
        if checker.check_ring(ring, final_product):
            final_product_scaffold = ring
            print(f"Final product contains {ring} scaffold")
            break

    # If no known tricyclic scaffold is found, check for bicyclic scaffolds that might be part of synthesis
    bicyclic_rings = [
        "indole",
        "quinoline",
        "isoquinoline",
        "naphthalene",
        "benzothiophene",
        "benzofuran",
        "indazole",
        "benzimidazole",
        "benzothiazole",
        "benzoxazole",
    ]

    if not final_product_scaffold:
        # Check for bicyclic scaffolds
        for ring in bicyclic_rings:
            if checker.check_ring(ring, final_product):
                final_product_scaffold = ring
                print(f"Final product contains bicyclic {ring} scaffold")
                break

    # If still no scaffold found, check if it has at least 2 connected rings
    if not final_product_scaffold:
        if ring_info.NumRings() >= 2:
            print("Final product has at least 2 rings, using general bicyclic/tricyclic detection")
            final_product_scaffold = "general_polycyclic"
        else:
            print("Final product does not have a suitable polycyclic scaffold")
            return False

    # Now trace the synthesis route to check if the scaffold is preserved
    scaffold_preserved = [True]  # Use a list to allow modification in nested function

    def check_scaffold_preservation(node, depth=0):
        if not scaffold_preserved[0]:
            return  # Early termination if we already know scaffold is not preserved

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check if the product has the scaffold
            product_has_scaffold = False
            if final_product_scaffold not in ["general_polycyclic"]:
                product_has_scaffold = checker.check_ring(final_product_scaffold, product)
            else:
                # For general polycyclic, check number of rings
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    product_ring_info = product_mol.GetRingInfo()
                    product_has_scaffold = product_ring_info.NumRings() >= 2

            # Check if at least one reactant has the scaffold
            reactant_has_scaffold = False
            for reactant in reactants:
                if final_product_scaffold not in ["general_polycyclic"]:
                    if checker.check_ring(final_product_scaffold, reactant):
                        reactant_has_scaffold = True
                        break
                else:
                    # For general polycyclic, check number of rings
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        reactant_ring_info = reactant_mol.GetRingInfo()
                        if reactant_ring_info.NumRings() >= 2:
                            reactant_has_scaffold = True
                            break

            # If product has scaffold but no reactant does, scaffold is not preserved
            if product_has_scaffold and not reactant_has_scaffold:
                print(
                    f"Scaffold not preserved at depth {depth}: Polycyclic scaffold created in this step"
                )
                print(f"Reaction: {rsmi}")
                scaffold_preserved[0] = False
                return

            # If product doesn't have scaffold but a reactant does, scaffold is not preserved
            if not product_has_scaffold and reactant_has_scaffold:
                print(
                    f"Scaffold not preserved at depth {depth}: Polycyclic scaffold broken in this step"
                )
                print(f"Reaction: {rsmi}")
                scaffold_preserved[0] = False
                return

        # Continue traversing the route
        for child in node.get("children", []):
            check_scaffold_preservation(child, depth + 1)

    check_scaffold_preservation(route)

    if scaffold_preserved[0]:
        print(f"Polycyclic scaffold ({final_product_scaffold}) preserved throughout synthesis")
    else:
        print(f"Polycyclic scaffold ({final_product_scaffold}) not preserved throughout synthesis")

    return scaffold_preserved[0]
