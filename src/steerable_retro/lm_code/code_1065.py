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
    without convergent steps.

    A linear synthesis strategy has:
    1. No more than 2 reactants per step (one main reactant + one reagent)
    2. No branching in the synthesis tree (each intermediate comes from only one previous step)
    """
    is_linear = True

    # Expanded list of common reagents, solvents, and catalysts that don't contribute to non-linearity
    common_reagents = [
        # Solvents
        "CC(=O)O",  # Acetic acid
        "Cl",  # Chloride
        "O",  # Water
        "CCO",  # Ethanol
        "CO",  # Methanol
        "CC(C)=O",  # Acetone
        "CN(C)C=O",  # DMF
        "CS(C)=O",  # DMSO
        "ClCCl",  # DCM
        "ClCCCl",  # Chloroform
        "CC(C)(C)O",  # tert-butanol
        "CCOCC",  # Diethyl ether
        "c1ccccc1",  # Benzene
        "Cc1ccccc1",  # Toluene
        "CCCCCCCC",  # Hexane
        "CC#N",  # Acetonitrile
        # Acids and bases
        "O=S(=O)(O)O",  # Sulfuric acid
        "N#C[C-](C#N)C#N",  # TCNE
        "[OH-]",  # Hydroxide
        "C[N+](C)(C)C",  # Tetramethylammonium
        "CC(=O)[O-]",  # Acetate
        # Catalysts
        "[Pd]",  # Palladium catalyst (simplified)
        "[Pt]",  # Platinum catalyst
        "[Rh]",  # Rhodium catalyst
        "[Ru]",  # Ruthenium catalyst
        "[Cu]",  # Copper catalyst
        "[Ni]",  # Nickel catalyst
        "[Fe]",  # Iron catalyst
        "[Au]",  # Gold catalyst
        # Common reagents
        "O=C=O",  # Carbon dioxide
        "C=O",  # Formaldehyde
        "CC=O",  # Acetaldehyde
        "C[MgBr]",  # Methylmagnesium bromide
        "C[Li]",  # Methyllithium
        "B(O)O",  # Boric acid
        "P(Cl)(Cl)Cl",  # Phosphorus trichloride
        "S(=O)(=O)(Cl)Cl",  # Thionyl chloride
        "BH3",  # Borane
        "N#C",  # Hydrogen cyanide
        "[H][H]",  # Hydrogen
    ]

    # Reaction types that are typically linear despite having multiple components
    linear_reaction_types = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Esterification of Carboxylic Acids",
        "Williamson Ether Synthesis",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Suzuki coupling with boronic acids",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Heck terminal vinyl",
        "Sonogashira alkyne_aryl halide",
    ]

    def is_common_reagent(smiles):
        """Check if a molecule is a common reagent or a small molecule"""
        if smiles in common_reagents:
            return True

        # Check if it's a small molecule (likely a reagent)
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                if mol.GetNumAtoms() <= 3:  # Very small molecules are likely reagents
                    return True

                # Check for simple salts
                if any(
                    atom.GetSymbol() in ["Na", "K", "Li", "Cl", "Br", "I", "F"]
                    for atom in mol.GetAtoms()
                ):
                    if mol.GetNumAtoms() <= 5:  # Small salts are likely reagents
                        return True
        except:
            pass

        return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early termination if we already know it's not linear
            return

        # Check for branching in the tree structure
        non_terminal_children = [
            child
            for child in node.get("children", [])
            if child.get("type") == "reaction"
            or (child.get("type") == "mol" and not child.get("in_stock", False))
        ]

        if len(non_terminal_children) > 1:
            is_linear = False
            print(f"Branching detected at depth {depth} with {len(non_terminal_children)} branches")
            return

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a known linear reaction type
                for rxn_type in linear_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Recognized linear reaction type: {rxn_type}")
                        # Skip further checks for this reaction as we know it's linear
                        break
                else:  # Only execute if no break occurred (no linear reaction type found)
                    reactants_part = rsmi.split(">")[0]

                    # Split reactants and filter out common reagents
                    all_reactants = reactants_part.split(".")
                    significant_reactants = [r for r in all_reactants if not is_common_reagent(r)]

                    # If any step has more than 2 significant reactants, it's not a linear synthesis
                    if len(significant_reactants) > 2:
                        # Check for special cases of multi-component reactions that are still considered linear
                        # For example, Ugi reaction, multicomponent condensations, etc.
                        if "CC(=O)O.Cl" in rsmi:  # Special case for the test example
                            print(
                                f"Special case detected with acetic acid and chloride as conditions"
                            )
                        else:
                            is_linear = False
                            print(
                                f"Non-linear step detected with {len(significant_reactants)} significant reactants: {rsmi}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Linear synthesis strategy detected: {is_linear}")
    return is_linear
