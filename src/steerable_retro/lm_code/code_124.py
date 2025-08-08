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


def main(route):
    """
    This function detects a convergent synthesis strategy where two complex fragments
    are joined in the final step.
    """
    is_convergent = False

    # Find all synthesized intermediates in the route
    synthesized_intermediates = set()

    def collect_intermediates(node):
        if node["type"] == "mol" and not node.get("in_stock", False):
            # Store both the original SMILES and a canonicalized version without atom mapping
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                synthesized_intermediates.add(node["smiles"])
                # Also add canonicalized SMILES without atom mapping
                canonical_smiles = Chem.MolToSmiles(mol)
                synthesized_intermediates.add(canonical_smiles)
        for child in node.get("children", []):
            collect_intermediates(child)

    collect_intermediates(route)
    print(f"Found {len(synthesized_intermediates)} synthesized intermediate entries in the route")

    # Find the final reaction (the first reaction node in retrosynthetic traversal)
    final_reaction = None

    def find_final_reaction(node):
        nonlocal final_reaction
        if node["type"] == "reaction" and final_reaction is None:
            final_reaction = node
        for child in node.get("children", []):
            find_final_reaction(child)

    find_final_reaction(route)

    if final_reaction and "rsmi" in final_reaction.get("metadata", {}):
        rsmi = final_reaction["metadata"]["rsmi"]
        print(f"Analyzing final reaction: {rsmi}")

        reactants = rsmi.split(">")[0].split(".")
        product = rsmi.split(">")[-1]

        print(f"Found {len(reactants)} reactants")

        # If we have at least 2 reactants in the final step
        if len(reactants) >= 2:
            # Check if both reactants are complex (not simple reagents)
            complex_synthesized_reactants = 0

            # Check each reactant in the final step
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    # Get canonicalized SMILES without atom mapping
                    canonical_reactant = Chem.MolToSmiles(mol)
                    # Also create a version with atom mapping removed
                    mol_no_mapping = Chem.MolFromSmiles(reactant)
                    for atom in mol_no_mapping.GetAtoms():
                        atom.SetAtomMapNum(0)
                    unmapped_reactant = Chem.MolToSmiles(mol_no_mapping)

                    # Define complexity: more than 7 atoms or contains a ring
                    atom_count = mol.GetNumAtoms()
                    ring_info = mol.GetRingInfo()
                    has_ring = ring_info.NumRings() > 0

                    print(f"Reactant: {reactant}")
                    print(f"Canonical: {canonical_reactant}")
                    print(f"Unmapped: {unmapped_reactant}")
                    print(f"Atoms: {atom_count}, Has ring: {has_ring}")

                    # Check if this reactant is a synthesized intermediate
                    # by checking all possible SMILES representations
                    is_intermediate = False

                    # Check if this reactant is a child of the final reaction
                    for child in final_reaction.get("children", []):
                        if child["type"] == "mol" and not child.get("in_stock", False):
                            child_mol = Chem.MolFromSmiles(child["smiles"])
                            if child_mol:
                                child_canonical = Chem.MolToSmiles(child_mol)
                                print(
                                    f"Comparing with child: {child['smiles']} (canonical: {child_canonical})"
                                )

                                # Try to match with the child's SMILES
                                if (
                                    reactant == child["smiles"]
                                    or canonical_reactant == child["smiles"]
                                    or unmapped_reactant == child["smiles"]
                                    or reactant == child_canonical
                                    or canonical_reactant == child_canonical
                                    or unmapped_reactant == child_canonical
                                ):
                                    is_intermediate = True
                                    print(f"Match found with child node!")
                                    break

                    # If not found in direct children, check all synthesized intermediates
                    if not is_intermediate:
                        if (
                            reactant in synthesized_intermediates
                            or canonical_reactant in synthesized_intermediates
                            or unmapped_reactant in synthesized_intermediates
                        ):
                            is_intermediate = True
                            print(f"Match found in synthesized intermediates!")

                    print(f"Is synthesized intermediate: {is_intermediate}")

                    # Count as complex synthesized reactant if it meets criteria
                    if (atom_count > 7 or has_ring) and is_intermediate:
                        complex_synthesized_reactants += 1
                        print(f"Found complex synthesized reactant: {reactant}")
                    # If it's complex but not found as synthesized, it might still be synthesized
                    # but with different SMILES representation - check children directly
                    elif atom_count > 7 or has_ring:
                        # For complex reactants, assume they're synthesized if they're not explicitly marked as in_stock
                        for child in final_reaction.get("children", []):
                            if child["type"] == "mol":
                                child_mol = Chem.MolFromSmiles(child["smiles"])
                                if child_mol and not child.get("in_stock", True):
                                    # Try structure matching instead of just SMILES comparison
                                    if mol.GetNumAtoms() == child_mol.GetNumAtoms():
                                        print(
                                            f"Complex reactant likely matches a synthesized intermediate"
                                        )
                                        complex_synthesized_reactants += 1
                                        print(f"Found complex synthesized reactant: {reactant}")
                                        break

            print(f"Complex synthesized reactants: {complex_synthesized_reactants}")
            if complex_synthesized_reactants >= 2:
                print("Convergent synthesis detected: final step joins two complex fragments")
                is_convergent = True

    return is_convergent
