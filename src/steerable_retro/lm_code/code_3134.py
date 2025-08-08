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
    This function detects if the synthetic route follows a linear synthesis strategy
    while preserving key structural motifs (triazolone, thioether, spiro system).
    """
    # Track when motifs first appear and if they're preserved
    triazolone_found = False
    thioether_found = False
    spiro_found = False

    # Track if motifs are preserved after first appearance
    triazolone_preserved = True
    thioether_preserved = True
    spiro_preserved = True

    # Track if synthesis is linear (no branching)
    is_linear = True
    reaction_count = 0

    # Collect all molecule nodes to check for branching
    all_mol_nodes = []

    def dfs_traverse(node, depth=0):
        nonlocal triazolone_found, thioether_found, spiro_found
        nonlocal triazolone_preserved, thioether_preserved, spiro_preserved
        nonlocal is_linear, reaction_count

        if node["type"] == "mol":
            # Add molecule node to list for branching check
            all_mol_nodes.append(node)

            # Check if this molecule node has multiple reactions (branching)
            if len(node.get("children", [])) > 1:
                print(f"Branching detected at molecule node, depth {depth}")
                is_linear = False

        elif node["type"] == "reaction":
            reaction_count += 1

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for triazolone pattern (triazole ring with carbonyl)
                has_triazolone = checker.check_ring("triazole", product) and checker.check_fg(
                    "Ketone", product
                )
                if has_triazolone:
                    print(f"Triazolone found at depth {depth}")
                    triazolone_found = True
                elif triazolone_found:
                    # Check if triazolone is preserved in reactants
                    reactants = rsmi.split(">")[0].split(".")
                    reactant_has_triazolone = any(
                        checker.check_ring("triazole", r) and checker.check_fg("Ketone", r)
                        for r in reactants
                    )
                    if not reactant_has_triazolone:
                        print(f"Triazolone pattern lost at depth {depth}")
                        triazolone_preserved = False

                # Check for thioether pattern
                has_thioether = checker.check_fg("Monosulfide", product)
                if has_thioether:
                    print(f"Thioether found at depth {depth}")
                    thioether_found = True
                elif thioether_found:
                    # Check if thioether is preserved in reactants
                    reactants = rsmi.split(">")[0].split(".")
                    reactant_has_thioether = any(
                        checker.check_fg("Monosulfide", r) for r in reactants
                    )
                    if not reactant_has_thioether:
                        print(f"Thioether pattern lost at depth {depth}")
                        thioether_preserved = False

                # Check for spiro system
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        ring_info = mol.GetRingInfo()
                        atom_rings = ring_info.AtomRings()
                        # Spiro system: two rings sharing exactly one atom
                        has_spiro = False
                        if len(atom_rings) >= 2:
                            for i in range(len(atom_rings)):
                                for j in range(i + 1, len(atom_rings)):
                                    shared_atoms = set(atom_rings[i]) & set(atom_rings[j])
                                    if len(shared_atoms) == 1:
                                        has_spiro = True
                                        break
                                if has_spiro:
                                    break

                        if has_spiro:
                            print(f"Spiro system found at depth {depth}")
                            spiro_found = True
                        elif spiro_found:
                            # Check if spiro is preserved in reactants
                            reactants = rsmi.split(">")[0].split(".")
                            reactant_has_spiro = False
                            for r in reactants:
                                r_mol = Chem.MolFromSmiles(r)
                                if r_mol:
                                    r_ring_info = r_mol.GetRingInfo()
                                    r_atom_rings = r_ring_info.AtomRings()
                                    if len(r_atom_rings) >= 2:
                                        for i in range(len(r_atom_rings)):
                                            for j in range(i + 1, len(r_atom_rings)):
                                                shared_atoms = set(r_atom_rings[i]) & set(
                                                    r_atom_rings[j]
                                                )
                                                if len(shared_atoms) == 1:
                                                    reactant_has_spiro = True
                                                    break
                                            if reactant_has_spiro:
                                                break
                                if reactant_has_spiro:
                                    break

                            if not reactant_has_spiro:
                                print(f"Spiro system lost at depth {depth}")
                                spiro_preserved = False
                except Exception as e:
                    print(f"Error analyzing spiro system: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Linear synthesis requires at least 2 reactions
    is_linear = is_linear and reaction_count >= 2

    # Check if all motifs were found and preserved
    motifs_preserved = (
        (not triazolone_found or triazolone_preserved)
        and (not thioether_found or thioether_preserved)
        and (not spiro_found or spiro_preserved)
    )

    result = is_linear and motifs_preserved

    print(f"Linear synthesis: {is_linear}")
    print(f"Triazolone found: {triazolone_found}, preserved: {triazolone_preserved}")
    print(f"Thioether found: {thioether_found}, preserved: {thioether_preserved}")
    print(f"Spiro found: {spiro_found}, preserved: {spiro_preserved}")
    print(f"Linear synthesis with preserved motifs: {result}")

    return result
