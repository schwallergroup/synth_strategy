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
    Detects if the synthetic route follows a convergent synthesis strategy
    with multiple fragment couplings rather than linear chain extension.
    """
    multi_fragment_reactions = 0
    debug = False  # Set to True to enable debug prints

    # Common reagents/catalysts/solvents to filter out as SMILES
    common_reagents = [
        "CCO",
        "C[O-]",
        "O=C(O)O",
        "[Pd]",
        "CC(=O)O",
        "CN(C)C",
        "CS(=O)(=O)O",
        "O",
        "N",
        "CC(C)C",
        "c1ccccc1",
        "Cc1ccccc1",
        "[Cs+]",
        "[K+]",
        "[Na+]",
        "B(O)O",
        "OB(O)O",
        "CC(=O)[O-]",
        "C1CCOC1",
        "ClCCl",
        "CN1CCCC1=O",
        "CCN(CC)CC",
        "CC#N",
        "C1CCCCC1",
        "O=S(=O)(O)O",
        "N#N",
        "C=O",
        "C(=O)O",
    ]

    # Convert common reagents to molecules for substructure matching
    common_reagent_mols = [
        Chem.MolFromSmiles(smiles) for smiles in common_reagents if Chem.MolFromSmiles(smiles)
    ]

    def is_common_reagent(smiles):
        """Check if a SMILES string matches common reagents/catalysts/solvents using RDKit"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Check if it's a very small molecule (likely a reagent)
        if mol.GetNumAtoms() < 3:
            return True

        # Check against our list of common reagents
        for reagent_mol in common_reagent_mols:
            if (
                mol.HasSubstructMatch(reagent_mol)
                and mol.GetNumAtoms() <= reagent_mol.GetNumAtoms() + 2
            ):
                return True

        return False

    def is_coupling_reaction(rsmi):
        """Check if a reaction is a coupling reaction"""
        # Check common coupling reactions
        coupling_reactions = [
            "Suzuki coupling with boronic acids",
            "Suzuki coupling with boronic esters",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
            "Stille reaction_aryl",
            "Stille reaction_vinyl",
            "Negishi coupling",
            "Sonogashira alkyne_aryl halide",
            "Sonogashira acetylene_aryl halide",
            "Heck terminal vinyl",
            "Hiyama-Denmark Coupling",
            "Kumada cross-coupling",
            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
            "Ullmann condensation",
        ]

        for reaction in coupling_reactions:
            if checker.check_reaction(reaction, rsmi):
                if debug:
                    print(f"Detected coupling reaction: {reaction}")
                return True

        return False

    def assess_fragment_complexity(mol):
        """Assess the complexity of a molecular fragment"""
        if mol is None:
            return 0

        # Basic complexity metrics
        atom_count = mol.GetNumAtoms()
        ring_count = rdMolDescriptors.CalcNumRings(mol)

        # Count functional groups as a measure of complexity
        fg_count = 0
        important_fgs = [
            "Aromatic halide",
            "Boronic acid",
            "Boronic ester",
            "Nitrile",
            "Carboxylic acid",
            "Ester",
            "Amide",
            "Amine",
            "Alcohol",
            "Ketone",
            "Aldehyde",
        ]

        for fg in important_fgs:
            if checker.check_fg(fg, Chem.MolToSmiles(mol)):
                fg_count += 1

        # Calculate complexity score
        complexity_score = 0

        # Size contribution
        if atom_count > 15:
            complexity_score += 2
        elif atom_count > 10:
            complexity_score += 1

        # Ring contribution
        complexity_score += min(ring_count, 3)

        # Functional group contribution
        complexity_score += min(fg_count, 2)

        if debug and complexity_score >= 2:
            print(
                f"Complex fragment: {Chem.MolToSmiles(mol)} - Score: {complexity_score} (Atoms: {atom_count}, Rings: {ring_count}, FGs: {fg_count})"
            )

        return complexity_score

    def dfs_traverse(node, depth=0):
        nonlocal multi_fragment_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Filter out common reagents
                significant_reactants = []
                for r in reactants:
                    if not is_common_reagent(r):
                        mol = Chem.MolFromSmiles(r)
                        if mol and mol.GetNumAtoms() > 3:  # Ensure it's not too small
                            significant_reactants.append(r)

                # Check if this is a coupling reaction
                is_coupling = is_coupling_reaction(rsmi)

                # Count complex fragments
                complex_reactants = 0
                complex_fragments = []

                for reactant in significant_reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    complexity = assess_fragment_complexity(mol)

                    if complexity >= 2:  # Threshold for considering a fragment complex
                        complex_reactants += 1
                        complex_fragments.append(reactant)

                # If a reaction has 2+ complex reactants or is a coupling with at least one complex reactant
                if complex_reactants >= 2 or (
                    is_coupling and complex_reactants >= 1 and len(significant_reactants) >= 2
                ):
                    multi_fragment_reactions += 1
                    if debug:
                        print(
                            f"Detected convergent synthesis step with {complex_reactants} complex fragments: {rsmi}"
                        )
                        print(f"Complex fragments: {complex_fragments}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Even one significant fragment coupling can indicate a convergent approach
    return multi_fragment_reactions >= 1
