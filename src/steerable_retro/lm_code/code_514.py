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
    This function detects a synthetic strategy involving formation of a tricyclic system.
    """
    forms_tricyclic_system = False

    def dfs_traverse(node, depth=0):
        nonlocal forms_tricyclic_system

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Get reactants as separate molecules
                reactants = reactants_part.split(".")

                # Check if product contains a tricyclic system
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Count rings in the product
                    ring_info = product_mol.GetRingInfo()
                    atom_rings = ring_info.AtomRings()
                    num_rings = len(atom_rings)

                    # Check if product has at least 3 rings
                    if num_rings >= 3:
                        # Check if rings are fused (form a tricyclic system)
                        fused_rings = False

                        # Check for shared atoms between rings (indicating fusion)
                        for i in range(num_rings):
                            for j in range(i + 1, num_rings):
                                # Convert ring atom tuples to sets for intersection
                                ring_i_atoms = set(atom_rings[i])
                                ring_j_atoms = set(atom_rings[j])

                                # If rings share atoms, they are fused
                                if len(ring_i_atoms.intersection(ring_j_atoms)) >= 2:
                                    fused_rings = True
                                    break
                            if fused_rings:
                                break

                        if fused_rings:
                            # Check if any reactant already has a tricyclic system
                            reactants_have_tricyclic = False
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    reactant_ring_info = reactant_mol.GetRingInfo()
                                    reactant_atom_rings = reactant_ring_info.AtomRings()
                                    reactant_num_rings = len(reactant_atom_rings)

                                    if reactant_num_rings >= 3:
                                        # Check if rings are fused in reactant
                                        for i in range(reactant_num_rings):
                                            for j in range(i + 1, reactant_num_rings):
                                                # Convert ring atom tuples to sets for intersection
                                                ring_i_atoms = set(reactant_atom_rings[i])
                                                ring_j_atoms = set(reactant_atom_rings[j])

                                                # If rings share atoms, they are fused
                                                if (
                                                    len(ring_i_atoms.intersection(ring_j_atoms))
                                                    >= 2
                                                ):
                                                    reactants_have_tricyclic = True
                                                    break
                                            if reactants_have_tricyclic:
                                                break

                                    if reactants_have_tricyclic:
                                        break

                            # If product has a tricyclic system but reactants don't, then it was formed
                            if not reactants_have_tricyclic:
                                print(f"Detected formation of a tricyclic system at depth {depth}")
                                forms_tricyclic_system = True

                            # Alternative approach: check for specific tricyclic ring systems
                            common_tricyclic_rings = [
                                "anthracene",
                                "phenanthrene",
                                "acridine",
                                "carbazole",
                                "dibenzofuran",
                                "dibenzothiophene",
                                "phenothiazine",
                                "phenoxazine",
                                "purine",
                                "xanthene",
                                "thioxanthene",
                                "pteridin",
                            ]

                            for ring_name in common_tricyclic_rings:
                                if checker.check_ring(ring_name, product):
                                    # Check if this ring was already in reactants
                                    ring_in_reactants = False
                                    for reactant in reactants:
                                        if checker.check_ring(ring_name, reactant):
                                            ring_in_reactants = True
                                            break

                                    if not ring_in_reactants:
                                        print(
                                            f"Detected formation of {ring_name} tricyclic system at depth {depth}"
                                        )
                                        forms_tricyclic_system = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return forms_tricyclic_system
