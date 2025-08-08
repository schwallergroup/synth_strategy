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
    This function detects if the synthesis follows a linear strategy without convergent steps.

    A linear synthesis strategy builds the target molecule through sequential steps where
    each reaction adds one main component to the growing molecule. In contrast, a convergent
    strategy involves combining multiple complex fragments in key steps.

    Returns:
        bool: True if the synthesis is linear, False if convergent steps are detected
    """
    is_linear = True

    # Common reagents that shouldn't count as complex fragments
    common_reagents = [
        "CC(=O)OC(=O)C",  # Acetic anhydride
        "CC(=O)O",  # Acetic acid
        "CC(=O)Cl",  # Acetyl chloride
        "CCO",  # Ethanol
        "CO",  # Methanol
        "C(=O)O",  # Formic acid
        "[OH-]",  # Hydroxide
        "O",  # Water
        "CC(=O)[O-]",  # Acetate
    ]

    def is_common_reagent(smiles):
        """Check if a molecule is a common reagent"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check if it's a small molecule (likely a reagent)
        if mol.GetNumHeavyAtoms() < 6:
            print(f"Small molecule detected as common reagent: {smiles}")
            return True

        # Check if it matches any common reagent patterns
        for reagent in common_reagents:
            reagent_mol = Chem.MolFromSmiles(reagent)
            if reagent_mol and mol.GetNumHeavyAtoms() <= reagent_mol.GetNumHeavyAtoms() + 2:
                print(f"Molecule matches common reagent pattern: {smiles}")
                return True

        return False

    def assess_complexity(smiles):
        """Assess the complexity of a molecule beyond just atom count"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0

        complexity = 0

        # Heavy atom count contributes to complexity
        heavy_atoms = mol.GetNumHeavyAtoms()
        complexity += heavy_atoms

        # Ring systems contribute significantly to complexity
        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumRings()
        complexity += ring_count * 3  # Rings add more complexity

        # Check for specific functional groups that indicate complexity
        if checker.check_fg("Aromatic alcohol", smiles):
            complexity += 2
            print(f"Found Aromatic alcohol in {smiles}")
        if checker.check_fg("Aromatic halide", smiles):
            complexity += 2
            print(f"Found Aromatic halide in {smiles}")
        if checker.check_fg("Ester", smiles):
            complexity += 2
            print(f"Found Ester in {smiles}")
        if checker.check_fg("Carboxylic acid", smiles):
            complexity += 2
            print(f"Found Carboxylic acid in {smiles}")
        if checker.check_fg("Aniline", smiles):
            complexity += 2
            print(f"Found Aniline in {smiles}")

        return complexity

    def is_coupling_reaction(rsmi):
        """Check if a reaction is a coupling reaction typically used in convergent synthesis"""
        try:
            if ">" not in rsmi or rsmi.count(">") < 2:
                print(f"Invalid reaction SMILES format: {rsmi}")
                return False

            is_suzuki = checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
            is_negishi = checker.check_reaction("Negishi coupling", rsmi)
            is_stille = checker.check_reaction("Stille reaction_aryl", rsmi)
            is_heck = checker.check_reaction("Heck terminal vinyl", rsmi)
            is_sonogashira = checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)

            result = is_suzuki or is_negishi or is_stille or is_heck or is_sonogashira

            if result:
                print(f"Coupling reaction detected: {rsmi}")

            return result
        except Exception as e:
            print(f"Error checking coupling reaction: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction" and is_linear:  # Skip if already determined to be non-linear
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    print(f"\nAnalyzing reaction at depth {depth}: {rsmi}")

                    if ">" not in rsmi or rsmi.count(">") < 2:
                        print(f"Invalid reaction SMILES format, skipping: {rsmi}")
                        return

                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")

                    # Check if this is a known coupling reaction
                    if is_coupling_reaction(rsmi):
                        significant_reactants = []
                        for r in reactants:
                            complexity = assess_complexity(r)
                            if not is_common_reagent(r) and complexity > 10:
                                significant_reactants.append((r, complexity))
                                print(
                                    f"Complex reactant in coupling reaction: {r} with complexity score {complexity}"
                                )

                        if len(significant_reactants) > 1:
                            # Check if one reactant is significantly more complex than others
                            significant_reactants.sort(key=lambda x: x[1], reverse=True)
                            primary_complexity = significant_reactants[0][1]
                            secondary_complexity = significant_reactants[1][1]

                            # If the most complex reactant is at least 2x more complex than the next one,
                            # it might still be considered linear
                            if primary_complexity < 2 * secondary_complexity:
                                print(f"Convergent coupling reaction detected at depth {depth}")
                                print(f"Primary reactant complexity: {primary_complexity}")
                                print(f"Secondary reactant complexity: {secondary_complexity}")
                                is_linear = False
                            else:
                                print(
                                    f"Coupling reaction at depth {depth} has a dominant reactant (complexity ratio: {primary_complexity/secondary_complexity:.2f})"
                                )
                    else:
                        # For non-coupling reactions, identify complex reactants
                        significant_reactants = []
                        for r in reactants:
                            if not is_common_reagent(r):
                                complexity = assess_complexity(r)
                                if complexity > 12:  # Higher threshold for non-coupling reactions
                                    significant_reactants.append((r, complexity))
                                    print(
                                        f"Complex reactant found: {r} with complexity score {complexity}"
                                    )

                        # If more than one complex reactant, check their relative complexity
                        if len(significant_reactants) > 1:
                            significant_reactants.sort(key=lambda x: x[1], reverse=True)
                            primary_complexity = significant_reactants[0][1]
                            secondary_complexity = significant_reactants[1][1]

                            # If the most complex reactant is at least 2.5x more complex than the next one,
                            # it might still be considered linear
                            if primary_complexity < 2.5 * secondary_complexity:
                                print(f"Convergent step detected at depth {depth}")
                                print(f"Primary reactant complexity: {primary_complexity}")
                                print(f"Secondary reactant complexity: {secondary_complexity}")
                                for idx, (smiles, complexity) in enumerate(significant_reactants):
                                    print(
                                        f"  Reactant {idx+1}: {smiles} (complexity score: {complexity})"
                                    )
                                is_linear = False
                            else:
                                print(
                                    f"Reaction at depth {depth} has a dominant reactant (complexity ratio: {primary_complexity/secondary_complexity:.2f})"
                                )
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"\nFinal determination: {'Linear' if is_linear else 'Convergent'} synthesis strategy")
    return is_linear
