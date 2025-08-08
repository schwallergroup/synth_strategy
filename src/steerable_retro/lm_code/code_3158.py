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
    Detects if the synthesis route includes an early-stage benzyl ether formation
    (formation of Ar-CH2-O-Ar linkage at depth 2+)
    """
    found_ether_formation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal found_ether_formation

        print(f"Traversing node at depth {current_depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction at depth {current_depth}: {rsmi}")

            # Check if this is at depth 2 or higher (early stage)
            if current_depth >= 2:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {current_depth}")
                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check if this is a Williamson Ether Synthesis or other ether formation reaction
                    is_williamson = checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    is_mitsunobu = checker.check_reaction("Mitsunobu aryl ether", rsmi)

                    # Check for benzyl group in reactants and ether in product
                    benzyl_halide_found = False
                    phenol_or_alcohol_found = False

                    for reactant in reactants:
                        print(f"Analyzing reactant: {reactant}")

                        # Check for benzyl halide (aromatic-CH2-X)
                        if "c" in reactant and (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or "Cl" in reactant
                        ):
                            # Check specifically for benzyl structure
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # More inclusive pattern for benzyl halides
                                benzyl_pattern = Chem.MolFromSmarts(
                                    "c-[CH2][F,Cl,Br,I,#9,#17,#35,#53]"
                                )
                                if reactant_mol.HasSubstructMatch(benzyl_pattern):
                                    benzyl_halide_found = True
                                    print(f"Found benzyl halide: {reactant}")
                                else:
                                    # Try a simpler pattern that might catch atom-mapped halides
                                    simple_benzyl_pattern = Chem.MolFromSmarts(
                                        "c-[CH2]~[#9,#17,#35,#53]"
                                    )
                                    if reactant_mol.HasSubstructMatch(simple_benzyl_pattern):
                                        benzyl_halide_found = True
                                        print(f"Found benzyl halide (simple pattern): {reactant}")

                        # Check for phenol or alcohol (for Mitsunobu or Williamson)
                        if (
                            checker.check_fg("Phenol", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or "OH" in reactant
                        ):
                            phenol_or_alcohol_found = True
                            print(f"Found phenol or alcohol: {reactant}")

                    # If we have the right reactants, check for ether in product
                    if benzyl_halide_found and phenol_or_alcohol_found:
                        print("Found both benzyl halide and phenol/alcohol reactants")

                        # Check if product contains ether
                        if checker.check_fg("Ether", product):
                            print("Product contains ether functional group")

                            # Verify benzyl ether formation in product
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                # Check for benzyl ether pattern
                                benzyl_ether_pattern = Chem.MolFromSmarts("c-[CH2]-[O]-c")
                                if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                                    print(
                                        f"Found early-stage benzyl ether formation at depth {current_depth}"
                                    )
                                    found_ether_formation = True
                                else:
                                    # Try a more general pattern
                                    general_ether_pattern = Chem.MolFromSmarts("c~[CH2]~[O]~c")
                                    if product_mol.HasSubstructMatch(general_ether_pattern):
                                        print(
                                            f"Found early-stage benzyl ether formation (general pattern) at depth {current_depth}"
                                        )
                                        found_ether_formation = True

                    # If we've already identified a Williamson or Mitsunobu reaction, do additional checks
                    if (is_williamson or is_mitsunobu) and not found_ether_formation:
                        print(
                            f"Found ether formation reaction: {'Williamson' if is_williamson else 'Mitsunobu'}"
                        )

                        # Check if product contains benzyl ether structure
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            benzyl_ether_pattern = Chem.MolFromSmarts("c~[CH2]~[O]~c")
                            if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                                print(
                                    f"Found early-stage benzyl ether formation in known reaction at depth {current_depth}"
                                )
                                found_ether_formation = True

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    print(f"Final result: {found_ether_formation}")
    return found_ether_formation
