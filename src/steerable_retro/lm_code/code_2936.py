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
    Detects if the synthetic route involves late-stage incorporation of a tetrahydropyran ring.
    """
    found_late_stage_thp = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_thp

        if node["type"] == "reaction":
            try:
                if "metadata" not in node or "rsmi" not in node["metadata"]:
                    print("Missing metadata or rsmi in reaction node")
                    for child in node.get("children", []):
                        dfs_traverse(child, depth + 1)
                    return

                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}")
                print(f"Reactants: {reactants}")
                print(f"Product: {product}")

                # Check if this is a late-stage reaction (depth ≤ 2)
                if depth <= 2:
                    # Check if product contains tetrahydropyran
                    if checker.check_ring("tetrahydropyran", product):
                        print("Product contains tetrahydropyran ring")

                        # Get atom indices of tetrahydropyran in product (with atom mapping)
                        product_thp_indices = checker.get_ring_atom_indices(
                            "tetrahydropyran", product
                        )
                        if not product_thp_indices:
                            print("Could not get tetrahydropyran atom indices in product")
                            for child in node.get("children", []):
                                dfs_traverse(child, depth + 1)
                            return

                        # Extract atom mapping numbers from product's tetrahydropyran
                        product_mol = Chem.MolFromSmiles(product)
                        if not product_mol:
                            print("Could not parse product SMILES")
                            for child in node.get("children", []):
                                dfs_traverse(child, depth + 1)
                            return

                        # Check if the tetrahydropyran in product is present in any reactant
                        thp_incorporated = False
                        thp_newly_formed = True

                        for reactant in reactants:
                            if checker.check_ring("tetrahydropyran", reactant):
                                print(f"Found tetrahydropyran in reactant: {reactant}")

                                # Get atom mapping numbers from reactant's tetrahydropyran
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if not reactant_mol:
                                    print(f"Could not parse reactant SMILES: {reactant}")
                                    continue

                                # Check if this is the same tetrahydropyran (by atom mapping)
                                # If all atoms in the tetrahydropyran ring have the same mapping in
                                # both product and reactant, then it's the same ring

                                # Extract atom mapping numbers from the tetrahydropyran in the reactant
                                reactant_thp_indices = checker.get_ring_atom_indices(
                                    "tetrahydropyran", reactant
                                )

                                # If we find the same tetrahydropyran structure in a reactant,
                                # it means the ring was not newly formed
                                if reactant_thp_indices:
                                    # Check if the atom mappings match between product and reactant
                                    # by comparing the atom mapping numbers

                                    # Get the atom mapping numbers from the product
                                    product_atom_maps = set()
                                    for ring in product_thp_indices:
                                        for atom_idx in ring:
                                            atom = product_mol.GetAtomWithIdx(atom_idx)
                                            if atom.HasProp("molAtomMapNumber"):
                                                map_num = atom.GetProp("molAtomMapNumber")
                                                product_atom_maps.add(int(map_num))

                                    # Get the atom mapping numbers from the reactant
                                    reactant_atom_maps = set()
                                    for ring in reactant_thp_indices:
                                        for atom_idx in ring:
                                            atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                            if atom.HasProp("molAtomMapNumber"):
                                                map_num = atom.GetProp("molAtomMapNumber")
                                                reactant_atom_maps.add(int(map_num))

                                    # If there's significant overlap in atom mapping numbers,
                                    # it's likely the same tetrahydropyran ring
                                    if (
                                        len(product_atom_maps.intersection(reactant_atom_maps)) >= 4
                                    ):  # At least 4 atoms match
                                        print(
                                            "Same tetrahydropyran ring found in reactant and product (by atom mapping)"
                                        )
                                        thp_newly_formed = False
                                        thp_incorporated = True
                                        break

                        # If tetrahydropyran is newly formed or incorporated in a late-stage reaction
                        if thp_newly_formed:
                            print(
                                "Found late-stage tetrahydropyran formation - ring is newly formed"
                            )
                            # Check for ring-forming reactions
                            if (
                                checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                                or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                                or checker.check_reaction("Diels-Alder", rsmi)
                            ):
                                print("Confirmed tetrahydropyran-forming reaction")
                                found_late_stage_thp = True
                            else:
                                print(
                                    "Tetrahydropyran appears but no specific ring-forming reaction detected"
                                )
                                # Still mark as found since the ring appears
                                found_late_stage_thp = True
                        elif thp_incorporated:
                            print(
                                "Found late-stage tetrahydropyran incorporation - existing ring is incorporated"
                            )
                            found_late_stage_thp = True
                        else:
                            print(
                                "Tetrahydropyran ring was already present in the molecule, not newly incorporated"
                            )

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_late_stage_thp
