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
    Detects if the synthesis follows a linear strategy (no convergent steps).

    A linear synthesis has at most one significant (non-in-stock or complex)
    reactant at each step. Complex molecules are defined as those with more
    than 7 heavy atoms that aren't marked as in_stock and aren't common reagents.
    """

    def is_common_reagent(smiles):
        """Check if a molecule is a common reagent that shouldn't count as a significant reactant"""
        # Common acids, bases, oxidants, reductants, etc.
        common_reagents = [
            "O=S(=O)(O)O",  # Sulfuric acid
            "O=C(O)C(F)(F)F",  # TFA
            "O=C(O)O",  # Carbonic acid
            "O=P(O)(O)O",  # Phosphoric acid
            "Cl[B-](Cl)(Cl)Cl.[K+]",  # Potassium tetrachloroborate
            "CC(=O)O",  # Acetic acid
            "O",  # Water
            "CO",  # Methanol
            "CCO",  # Ethanol
            "CC(C)O",  # Isopropanol
            "CC(O)C",  # Propanol
            "ClCCl",  # Dichloromethane
            "ClC(Cl)Cl",  # Chloroform
            "CC(=O)OC(=O)C",  # Acetic anhydride
            "O=S(=O)(O)S(=O)(=O)O",  # Sulfonic acid
            "N#CC#N",  # Malononitrile
            "C1CCOC1",  # THF
            "C1CCOCC1",  # THP
            "CN(C)C=O",  # DMF
            "CS(=O)C",  # DMSO
            "CC#N",  # Acetonitrile
            "c1ccccc1",  # Benzene
            "Cc1ccccc1",  # Toluene
            "CC(C)(C)O",  # tert-butanol
            "O=C=O",  # Carbon dioxide
            "N",  # Ammonia
            "Cl",  # Chlorine
            "Br",  # Bromine
            "I",  # Iodine
            "F",  # Fluorine
            "O=O",  # Oxygen
            "[H][H]",  # Hydrogen
        ]

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check if molecule is a common reagent by direct SMILES comparison
        for reagent in common_reagents:
            if Chem.MolToSmiles(Chem.MolFromSmiles(reagent)) == Chem.MolToSmiles(mol):
                return True

        # Check if molecule is a simple acid, base, or salt
        if mol.GetNumHeavyAtoms() <= 7:
            # Check for common functional groups that indicate reagents
            if (
                checker.check_fg("Carboxylic acid", smiles)
                or checker.check_fg("Sulfonic acid", smiles)
                or checker.check_fg("Primary amine", smiles)
                or checker.check_fg("Secondary amine", smiles)
                or checker.check_fg("Tertiary amine", smiles)
            ):
                return True

        return False

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            print(f"Examining reaction at depth {depth}")

            # Get molecule children (reactants)
            mol_children = [child for child in node.get("children", []) if child["type"] == "mol"]

            # Get reaction SMILES if available
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                print(f"  Reaction SMILES: {rsmi}")
                reactants_smiles = rsmi.split(">")[0].split(".")
            else:
                reactants_smiles = []
                print("  No reaction SMILES available")

            # Count significant reactants (not in_stock, complex, and not common reagents)
            significant_reactants = []

            for child in mol_children:
                child_smiles = child["smiles"]

                # Skip in-stock molecules (starting materials)
                if child.get("in_stock", False):
                    print(f"  Skipping in-stock reactant: {child_smiles}")
                    continue

                # Skip if it's a common reagent
                if is_common_reagent(child_smiles):
                    print(f"  Skipping common reagent: {child_smiles}")
                    continue

                # Check molecule complexity by atom count
                mol = Chem.MolFromSmiles(child_smiles)
                if mol and mol.GetNumHeavyAtoms() > 7:
                    # Verify this molecule is actually a reactant in the reaction SMILES
                    is_reactant = True
                    if reactants_smiles:
                        # Simple check - see if the molecule appears in the reactants
                        # This is a simplification; ideally we'd use atom mapping
                        is_reactant = any(
                            Chem.MolFromSmiles(r)
                            and Chem.MolToSmiles(Chem.MolFromSmiles(r)) == Chem.MolToSmiles(mol)
                            for r in reactants_smiles
                        )

                    if is_reactant:
                        significant_reactants.append(child_smiles)
                        print(
                            f"  Found significant reactant: {child_smiles} with {mol.GetNumHeavyAtoms()} heavy atoms"
                        )

            if len(significant_reactants) > 1:
                print(f"Convergent step detected at depth {depth} - not a linear synthesis")
                print(f"Reaction has {len(significant_reactants)} significant reactants:")
                for i, smiles in enumerate(significant_reactants):
                    print(f"  Reactant {i+1}: {smiles}")

                return False
            else:
                print(
                    f"  Linear step confirmed at depth {depth} with {len(significant_reactants)} significant reactants"
                )

        # Continue traversal
        for child in node.get("children", []):
            if not dfs_traverse(child, depth + 1):
                return False

        return True

    result = dfs_traverse(route)
    print(f"Final result: {'Linear synthesis' if result else 'Convergent synthesis'}")
    return result
