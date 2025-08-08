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
    (as opposed to convergent synthesis).

    A linear synthesis adds one component at a time to a growing molecule.
    A convergent synthesis combines complex fragments that were synthesized separately.
    """
    is_linear = True

    # Common reagents and catalysts that shouldn't count toward convergence
    common_reagents = [
        "O",
        "O=O",
        "[OH-]",
        "[H+]",
        "[Na+]",
        "[K+]",
        "Cl",
        "Br",
        "I",
        "CC(=O)O",
        "CCO",
        "CO",
        "CS(=O)(=O)O",
        "CN",
        "C1CCOC1",
        "ClC(Cl)(Cl)Cl",
        "CC(C)O",
        "O=S(=O)(O)O",
        "N#N",
        "N",
        "O=C=O",
        "C(=O)O",
        "C=O",
        "S",
        "P",
        "B",
        "Si",
        "CCN",
        "CC(C)N",
        "CN(C)C",
        "CC(N)=O",
        "CC#N",
    ]

    def is_reagent(smiles):
        """Check if a molecule is likely a reagent rather than a main reactant"""
        # Check against common reagents list
        if smiles in common_reagents:
            return True

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Small molecules are often reagents
        if mol.GetNumHeavyAtoms() <= 5:
            return True

        # Check for common reagent functional groups
        if (
            checker.check_fg("Triflate", smiles)
            or checker.check_fg("Tosylate", smiles)
            or checker.check_fg("Mesylate", smiles)
            or checker.check_fg("Acyl halide", smiles)
        ):
            return True

        return False

    def count_significant_atoms(mol):
        """Count atoms that are significant for determining structure"""
        count = 0
        has_ring = False

        # Check if molecule contains rings
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            has_ring = True

        for atom in mol.GetAtoms():
            # Skip hydrogen and common leaving groups
            if atom.GetSymbol() not in ["H", "F", "Cl", "Br", "I"]:
                # Atoms in rings are more significant
                if has_ring and atom.IsInRing():
                    count += 2  # Count ring atoms as more significant
                else:
                    count += 1

        return count

    def dfs_traverse(node):
        nonlocal is_linear

        if not is_linear:  # Early return if we already know it's not linear
            return

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants_smiles = reactants_part.split(".")

                # Get product molecule
                product_mol = Chem.MolFromSmiles(product_part)
                if not product_mol:
                    print(f"Warning: Could not parse product SMILES: {product_part}")
                    return

                product_size = product_mol.GetNumHeavyAtoms()

                # Filter out reagents and get main reactants
                main_reactants = []
                for reactant in reactants_smiles:
                    if not is_reagent(reactant):
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumHeavyAtoms() > 5:  # Increased threshold
                            main_reactants.append(reactant)

                # If more than one main reactant, check if it's truly convergent
                if len(main_reactants) > 1:
                    # Get atom maps from reactants
                    atom_contributions = {}
                    largest_contribution = 0

                    for reactant in main_reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        significant_atoms = count_significant_atoms(mol)

                        # Count mapped atoms
                        mapped_atoms = 0
                        for atom in mol.GetAtoms():
                            if atom.GetAtomMapNum() > 0:
                                mapped_atoms += 1

                        # If this reactant contributes significant structure, record it
                        if mapped_atoms >= 5 and significant_atoms >= 7:  # Increased thresholds
                            atom_contributions[reactant] = mapped_atoms
                            largest_contribution = max(largest_contribution, mapped_atoms)

                    # Check if multiple reactants contribute significant structure
                    if len(atom_contributions) > 1:
                        # Check if all contributions are significant relative to product size
                        significant_contributions = 0
                        for reactant, count in atom_contributions.items():
                            # A reactant is significant if it contributes at least 20% of product atoms
                            # and is not just a small addition to the largest fragment
                            if count >= product_size * 0.2 and count >= largest_contribution * 0.4:
                                significant_contributions += 1

                        # If multiple significant contributions, it's convergent
                        if significant_contributions > 1:
                            print(f"Convergent step detected in reaction: {rsmi}")
                            print(f"Multiple significant reactants: {atom_contributions}")
                            is_linear = False
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return is_linear
