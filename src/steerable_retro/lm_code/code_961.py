#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis follows a linear strategy
    (as opposed to convergent) by checking if each reaction has only
    one non-reagent reactant.
    """
    is_linear = True

    def is_common_reagent(smiles):
        """
        Identifies common reagents that shouldn't count toward convergence.
        """
        # Create RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Common reagents with small molecular weight
        if Descriptors.MolWt(mol) < 50:
            print(f"Reagent identified by small MW: {smiles}")
            return True

        # Check for common reagent functional groups
        common_reagent_groups = [
            "Triflate",
            "Mesylate",
            "Tosylate",
            "Magnesium halide",
            "Zinc halide",
            "Alkyl lithium",
            "Tin",
            "Silyl protective group",
            "TMS ether protective group",
            "Silane",
        ]

        for fg in common_reagent_groups:
            if checker.check_fg(fg, smiles):
                print(f"Reagent identified by functional group {fg}: {smiles}")
                return True

        # Check for phosphonium compounds (Wittig reagents)
        if "[P+]" in smiles or "P(C)(C)(C)(C)" in smiles or "P(=O)" in smiles:
            print(f"Phosphonium/phosphorus reagent identified: {smiles}")
            return True

        # Common reagent molecules - expanded list
        common_reagents = [
            "CCO",
            "CO",
            "O",
            "[OH-]",
            "[H+]",
            "Cl",
            "Br",
            "I",
            "F",
            "CC(=O)O",
            "CC(=O)[O-]",
            "CN",
            "CS(=O)(=O)O",
            "CS(=O)(=O)[O-]",
            "C[N+](C)(C)C",
            "B(O)O",
            "B(OH)3",
            "P(OCC)(OCC)OCC",
            "P(=O)(OCC)(OCC)OCC",
            "CC(C)O",
            "CCOC(=O)C",
            "CC(=O)C",
            "O=C=O",
            "N#N",
            "N=[N+]=[N-]",
            "[Na+]",
            "[K+]",
            "[Li+]",
            "c1ccccc1",
            "CC(=O)OC(=O)C",
            "ClC(=O)C",
            "BrC(=O)C",
            "IC(=O)C",
            "FC(=O)C",
            "O=S(=O)(O)O",
            "[O-]S(=O)(=O)[O-]",
            "N",
            "NN",
        ]

        # Normalize SMILES for comparison
        norm_smiles = Chem.MolToSmiles(mol)
        for reagent in common_reagents:
            reagent_mol = Chem.MolFromSmiles(reagent)
            if reagent_mol and Chem.MolToSmiles(reagent_mol) == norm_smiles:
                print(f"Common reagent molecule identified: {smiles}")
                return True

        # Check for low complexity molecules
        if Descriptors.NumRotatableBonds(mol) <= 2 and mol.GetNumAtoms() < 8:
            print(f"Simple molecule identified as reagent: {smiles}")
            return True

        return False

    def is_inherently_convergent_reaction(rsmi):
        """
        Checks if the reaction type is inherently convergent.
        """
        convergent_reaction_types = [
            "Suzuki coupling with boronic acids",
            "Suzuki coupling with boronic esters",
            "Negishi coupling",
            "Stille reaction_aryl",
            "Heck terminal vinyl",
            "Sonogashira alkyne_aryl halide",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
            "Ugi reaction",
            "A3 coupling",
        ]

        for rxn_type in convergent_reaction_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Inherently convergent reaction detected: {rxn_type}")
                return True

        return False

    def is_inherently_linear_reaction(rsmi):
        """
        Checks if the reaction type is inherently linear despite having multiple reactants.
        """
        linear_reaction_types = [
            "Wittig reaction",
            "Wittig",
            "Wittig reaction with triphenylphosphorane",
            "Wittig with Phosphonium",
            "Julia Olefination",
            "Grignard with CO2 to carboxylic acid",
            "Grignard from aldehyde to alcohol",
            "Grignard from ketone to alcohol",
        ]

        for rxn_type in linear_reaction_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Inherently linear reaction detected despite multiple reactants: {rxn_type}")
                return True

        return False

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and is_linear:  # Skip if already non-linear
            # Extract reactants
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Check if reaction type is inherently convergent
                if is_inherently_convergent_reaction(rsmi):
                    print(f"Non-linear step detected: inherently convergent reaction")
                    is_linear = False
                    return

                # Check if reaction type is inherently linear despite multiple reactants
                if is_inherently_linear_reaction(rsmi):
                    print(f"Linear step maintained: reaction is inherently linear")
                    # Continue with traversal
                else:
                    # Count significant reactants (excluding common reagents)
                    significant_reactants = []
                    for reactant in reactants:
                        if not is_common_reagent(reactant):
                            significant_reactants.append(reactant)

                    if len(significant_reactants) > 1:
                        print(
                            f"Non-linear step detected with {len(significant_reactants)} significant reactants"
                        )
                        for r in significant_reactants:
                            print(f"  - {r}")
                        is_linear = False
                        return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return is_linear
