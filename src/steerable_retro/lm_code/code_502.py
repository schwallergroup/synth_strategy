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
    This function detects a linear synthesis strategy where each step builds on a single precursor
    (as opposed to convergent synthesis where multiple fragments are combined).
    """
    # Track if all reactions have only one product
    all_steps_linear = True
    reaction_count = 0

    # List of reaction types that are considered linear despite having multiple reactants
    linear_reaction_types = [
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Esterification of Carboxylic Acids",
        "Williamson Ether Synthesis",
        "S-alkylation of thiols",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Sonogashira alkyne_aryl halide",
        "Heck terminal vinyl",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Alkylation of amines",
    ]

    def dfs_traverse(node):
        nonlocal all_steps_linear, reaction_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reaction_count += 1

            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            products = products_part.split(".")

            # If there are multiple reactants and only one product, check if it's a known linear reaction type
            if len(reactants) > 1 and len(products) == 1:
                is_known_linear_reaction = False

                # Check if this is a known linear reaction type
                for reaction_type in linear_reaction_types:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected linear reaction type: {reaction_type}")
                        is_known_linear_reaction = True
                        break

                # If not identified by reaction type, check for specific functional group transformations
                if not is_known_linear_reaction:
                    # Check for N-alkylation by examining the reaction
                    if checker.check_reaction("Alkylation of amines", rsmi):
                        print(f"Detected amine alkylation (linear reaction)")
                        is_known_linear_reaction = True

                    # Check for N-arylation
                    elif checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    ):
                        print(f"Detected N-arylation (linear reaction)")
                        is_known_linear_reaction = True

                    # Check for sulfonamide formation
                    elif checker.check_fg("Sulfonamide", products[0]):
                        print(f"Detected sulfonamide formation (linear reaction)")
                        is_known_linear_reaction = True

                    # Check for protecting group additions
                    elif any(
                        checker.check_fg(fg, products[0])
                        for fg in ["Boc", "TMS ether protective group", "Silyl protective group"]
                    ):
                        print(f"Detected protecting group addition (linear reaction)")
                        is_known_linear_reaction = True

                    # Check for heterocyclic N-alkylation
                    else:
                        # Try to identify N-alkylation by examining the product
                        product_mol = Chem.MolFromSmiles(products[0])
                        if product_mol:
                            # Check if the product contains heterocycles
                            has_heterocycle = False
                            for atom in product_mol.GetAtoms():
                                if atom.GetSymbol() == "N" and atom.IsInRing():
                                    has_heterocycle = True
                                    break

                            if has_heterocycle:
                                # Check if any reactant contains an alkyl halide or similar leaving group
                                for r in reactants:
                                    if (
                                        checker.check_fg("Primary halide", r)
                                        or checker.check_fg("Secondary halide", r)
                                        or checker.check_fg("Tertiary halide", r)
                                    ):
                                        print(
                                            f"Detected heterocyclic N-alkylation (linear reaction)"
                                        )
                                        is_known_linear_reaction = True
                                        break

                    # Check for metal-catalyzed reactions which are often linear despite multiple reactants
                    if not is_known_linear_reaction:
                        # More robust catalyst detection
                        catalyst_patterns = [
                            "Pd",
                            "Pt",
                            "Rh",
                            "Ru",
                            "Cu",
                            "Ni",
                            "Fe",
                            "Co",
                            "Ir",
                            "Au",
                            "Ag",
                        ]

                        # Check if any reactant contains catalyst-like patterns
                        catalyst_present = False
                        for r in reactants:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol:
                                # Check for metal atoms in the molecule
                                for atom in r_mol.GetAtoms():
                                    if atom.GetSymbol() in catalyst_patterns:
                                        catalyst_present = True
                                        break

                        if catalyst_present:
                            print(f"Detected metal-catalyzed reaction (likely linear)")
                            is_known_linear_reaction = True

                # If not a known linear reaction, it's potentially convergent
                if not is_known_linear_reaction:
                    # One last check - see if one reactant is significantly larger than the others
                    # This might indicate a main substrate with reagents
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]
                    if reactant_mols:
                        atom_counts = [mol.GetNumAtoms() for mol in reactant_mols]
                        max_atoms = max(atom_counts)
                        # If one reactant has significantly more atoms than others, it's likely the main substrate
                        if max_atoms > 10 and max_atoms > 2 * sum(atom_counts) / len(atom_counts):
                            print(f"Detected main substrate with reagents (likely linear)")
                            is_known_linear_reaction = True

                # If still not identified as linear, mark as convergent
                if not is_known_linear_reaction:
                    print(f"Detected potentially convergent step: {rsmi}")
                    all_steps_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Only consider it a linear strategy if there are multiple reactions and all are linear
    strategy_present = reaction_count > 1 and all_steps_linear
    print(f"Linear synthesis strategy detected: {strategy_present}")
    return strategy_present
