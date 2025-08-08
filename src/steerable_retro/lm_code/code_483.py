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
    Detects if the synthesis maintains a specific heterocyclic core (thiophene-fused pyridine)
    throughout the synthesis.
    """
    # Track if the core is present in final product and intermediates
    core_present_in_final = False
    core_present_in_intermediates = False

    # Track the path of molecules that contain the core
    core_path = []

    def has_thiophene_fused_pyridine(smiles):
        """Check if a molecule contains a thiophene-fused pyridine core."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                print(f"Could not parse SMILES: {smiles}")
                return False

            # Direct check for thiophene-fused pyridine (thienopyridine) using SMARTS patterns
            # Pattern 1: Thieno[2,3-b]pyridine
            pattern1 = Chem.MolFromSmarts("c1csc2c1ncccc2")
            # Pattern 2: Thieno[3,2-b]pyridine
            pattern2 = Chem.MolFromSmarts("c1ccnc2c1scc2")
            # Pattern 3: Thieno[2,3-c]pyridine
            pattern3 = Chem.MolFromSmarts("c1cncc2c1scc2")
            # Pattern 4: Thieno[3,2-c]pyridine
            pattern4 = Chem.MolFromSmarts("c1cncc2c1csc2")

            if (
                mol.HasSubstructMatch(pattern1)
                or mol.HasSubstructMatch(pattern2)
                or mol.HasSubstructMatch(pattern3)
                or mol.HasSubstructMatch(pattern4)
            ):
                print(f"Found thiophene-fused pyridine in: {smiles}")
                return True

            # Fallback method: check if both rings are present and fused
            has_thiophene = checker.check_ring("thiophene", smiles)
            has_pyridine = checker.check_ring("pyridine", smiles)

            if not (has_thiophene and has_pyridine):
                return False

            # Check for fusion using RDKit's ring info
            ring_info = mol.GetRingInfo()
            ring_atoms = ring_info.AtomRings()

            # Convert to sets for intersection checking
            ring_atom_sets = [set(ring) for ring in ring_atoms]

            # Check if any two rings share at least 2 atoms (fusion)
            for i in range(len(ring_atom_sets)):
                for j in range(i + 1, len(ring_atom_sets)):
                    if len(ring_atom_sets[i].intersection(ring_atom_sets[j])) >= 2:
                        # Found fused rings, now check if they are thiophene and pyridine
                        print(f"Found fused rings in: {smiles}")
                        return True

            return False
        except Exception as e:
            print(f"Error checking for thiophene-fused pyridine: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal core_present_in_final, core_present_in_intermediates

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            has_core = has_thiophene_fused_pyridine(smiles)

            if has_core:
                # Add to core path with depth information
                core_path.append((smiles, depth))
                print(f"Added to core path: {smiles} at depth {depth}")

                if depth == 0:
                    core_present_in_final = True
                    print("Core present in final product")
                else:
                    core_present_in_intermediates = True
                    print(f"Core present in intermediate at depth {depth}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # For reaction nodes, verify the core is maintained through the reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has the core
                product_has_core = has_thiophene_fused_pyridine(product)

                # Check if at least one reactant has the core
                reactant_has_core = any(has_thiophene_fused_pyridine(r) for r in reactants)

                if product_has_core and not reactant_has_core:
                    # If product has core but no reactant does, core was created in this step
                    print(f"Core was created in reaction: {rsmi}")
                elif product_has_core and reactant_has_core:
                    # Core is maintained through this reaction
                    print(f"Core is maintained through reaction: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a continuous path of core-containing molecules
    if len(core_path) >= 2:
        # Sort by depth (ascending)
        core_path.sort(key=lambda x: x[1])
        print(f"Core path found: {core_path}")

        # Check if the path includes both the final product and at least one intermediate
        depths = [d for _, d in core_path]
        if 0 in depths and any(d > 0 for d in depths):
            print("Core maintained throughout synthesis")
            return True

    result = core_present_in_final and core_present_in_intermediates
    print(f"Final result: {result}")
    return result
