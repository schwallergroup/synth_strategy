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
    This function detects if the synthesis follows a linear strategy (no convergent steps).
    A linear synthesis has no steps where multiple significant building blocks are combined.
    """
    is_linear = True

    # Common reagents and small molecules that shouldn't count toward convergence
    common_reagents = [
        "O",
        "CO",
        "CCO",
        "C(=O)O",
        "CC(=O)O",
        "CN",
        "C[NH]",
        "C[NH2]",
        "CS(=O)(=O)O",
        "CS(=O)(=O)[O-]",
        "C[N+](C)(C)C",
        "C[NH3+]",
        "Cl",
        "Br",
        "I",
        "F",
        "[OH-]",
        "[H+]",
        "[Na+]",
        "[K+]",
        "[Cl-]",
        "CC(=O)[O-]",
        "CC(=O)O[Na]",
        "CC(=O)O[K]",
        "C[O-]",
        "C[O-][Na+]",
        "C[O-][K+]",
        "O=C=O",
        "N#N",
        "[H][H]",
        "O=O",
        "N#[C-]",
        "[BH4-]",
        "B(O)O",
        "B(O)(O)O",
        "P(O)(O)O",
        "S(=O)(=O)(O)O",
        "CC(=O)OC(=O)C",
        "C[S](=O)(=O)[O-]",
        "C[S](=O)(=O)[O-][Na+]",
        "C[S](=O)(=O)[O-][K+]",
        "C[S](=O)(=O)Cl",
        "C[S](=O)(=O)F",
    ]

    def is_reagent(smiles):
        """Check if a molecule is likely a reagent rather than a building block"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"  Warning: Could not parse SMILES: {smiles}")
                return True

            # Small molecules are likely reagents
            if mol.GetNumHeavyAtoms() <= 6:
                return True

            # Check if it's a common reagent
            if smiles in common_reagents:
                return True

            # Check for common functional groups that indicate reagents
            if (
                checker.check_fg("Triflate", smiles)
                or checker.check_fg("Tosylate", smiles)
                or checker.check_fg("Mesylate", smiles)
                or checker.check_fg("Magnesium halide", smiles)
                or checker.check_fg("Zinc halide", smiles)
                or checker.check_fg("Tin", smiles)
                or checker.check_fg("Silane", smiles)
                or checker.check_fg("Alkyl lithium", smiles)
                or checker.check_fg("Aryl lithium", smiles)
                or checker.check_fg("Boronic acid", smiles)
                or checker.check_fg("Boronic ester", smiles)
                or checker.check_fg("Phosphate ester", smiles)
            ):
                return True

            return False
        except Exception as e:
            print(f"  Error in is_reagent: {e} for SMILES: {smiles}")
            return True  # If we can't process it, assume it's a reagent

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early return if we already know it's not linear
            return

        if node["type"] == "reaction":
            try:
                # Extract reactants
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Count significant building blocks (non-reagents)
                significant_reactants = []
                for r_smi in reactants_smiles:
                    if not is_reagent(r_smi):
                        significant_reactants.append(r_smi)
                        print(f"  Significant reactant: {r_smi}")
                    else:
                        print(f"  Reagent (not counted): {r_smi}")

                if len(significant_reactants) > 1:
                    is_linear = False
                    print(
                        f"Convergent step detected with {len(significant_reactants)} significant reactants"
                    )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy confirmed")
    else:
        print("Convergent synthesis strategy detected")

    return is_linear
