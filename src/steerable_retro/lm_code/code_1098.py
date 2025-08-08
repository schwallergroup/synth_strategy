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
    Detects if the synthesis follows a convergent approach where two complex fragments
    are joined together in a late-stage reaction.
    """
    convergent_detected = False

    # Common reagents and solvents to filter out
    common_reagents = [
        "O",
        "[OH-]",
        "[OH2+]",
        "O=O",
        "[O-]",
        "[O]",
        "[OH3+]",
        "[Na+]",
        "[K+]",
        "[Li+]",
        "[Mg+2]",
        "[Zn+2]",
        "[Cu+]",
        "[Cu+2]",
        "C1CCOC1",
        "CN",
        "CCO",
        "CO",
        "CC(=O)O",
        "CC(C)=O",
        "CS(=O)(=O)O",
        "ClCCl",
        "BrCCBr",
        "c1ccccc1",
        "N#N",
        "[Pd]",
        "[Pt]",
        "[Rh]",
        "C(=O)O",
        "C=O",
        "C#N",
        "C(C)(C)O",
        "C(F)(F)F",
        "C(Cl)(Cl)Cl",
    ]

    def is_common_reagent(smiles):
        """Check if a molecule is a common reagent or solvent"""
        # Direct match with common reagents list
        if smiles in common_reagents:
            return True

        # Check for very small molecules (likely reagents)
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumAtoms() <= 3:
            return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_detected

        if node["type"] == "reaction" and depth <= 3:  # Focus on late-stage reactions
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Filter out common reagents
                significant_reactants = [r for r in reactants if not is_common_reagent(r)]

                # Check if we have multiple significant reactants
                if len(significant_reactants) >= 2:
                    # Check if this is a coupling reaction
                    is_coupling = any(
                        [
                            checker.check_reaction("Suzuki coupling with boronic acids", rsmi),
                            checker.check_reaction("Suzuki coupling with boronic esters", rsmi),
                            checker.check_reaction("Negishi coupling", rsmi),
                            checker.check_reaction("Stille reaction_aryl", rsmi),
                            checker.check_reaction("Heck terminal vinyl", rsmi),
                            checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi),
                            checker.check_reaction(
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                            ),
                            checker.check_reaction("Ullmann condensation", rsmi),
                            checker.check_reaction(
                                "Wittig reaction with triphenylphosphorane", rsmi
                            ),
                            checker.check_reaction("Mitsunobu esterification", rsmi),
                            checker.check_reaction("Amidation", rsmi),
                            checker.check_reaction("Esterification of Carboxylic Acids", rsmi),
                            checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                rsmi,
                            ),
                        ]
                    )

                    # Count complex reactants
                    complex_reactants = []
                    for reactant in significant_reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            atom_count = mol.GetNumAtoms()
                            ring_count = mol.GetRingInfo().NumRings()
                            print(
                                f"  Reactant: {reactant}, Atoms: {atom_count}, Rings: {ring_count}"
                            )

                            # Define complexity as having significant structure
                            if atom_count >= 7 or ring_count >= 1:
                                complex_reactants.append(reactant)

                    print(
                        f"  Complex reactants: {len(complex_reactants)}, Is coupling: {is_coupling}"
                    )

                    # Check for amide formation, esterification, or other fragment-joining reactions
                    if not is_coupling:
                        is_joining_reaction = any(
                            [
                                checker.check_reaction(
                                    "Carboxylic acid with primary amine to amide", rsmi
                                ),
                                checker.check_reaction("Esterification of Carboxylic Acids", rsmi),
                                checker.check_reaction(
                                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                                ),
                                checker.check_reaction("Williamson Ether Synthesis", rsmi),
                            ]
                        )
                        if is_joining_reaction:
                            is_coupling = True

                    # Convergent synthesis detected if:
                    # 1. We have at least 2 complex reactants, OR
                    # 2. It's a coupling/joining reaction with at least 2 reactants where one is complex
                    if len(complex_reactants) >= 2 or (
                        is_coupling
                        and len(complex_reactants) >= 1
                        and len(significant_reactants) >= 2
                    ):
                        convergent_detected = True
                        print(
                            f"Convergent synthesis detected at depth {depth} with {len(complex_reactants)} complex fragments"
                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_detected
