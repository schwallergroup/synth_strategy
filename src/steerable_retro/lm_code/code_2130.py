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
    (each reaction has only one non-reagent reactant).

    Linear synthesis strategies typically involve:
    1. One main reactant and building blocks/reagents in each step
    2. Common coupling reactions like Suzuki, Sonogashira, Buchwald-Hartwig are considered linear
       even though they have two significant reactants
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a common coupling reaction (considered linear by convention)
            is_coupling_reaction = (
                # Suzuki couplings
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or
                # Sonogashira couplings
                checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
                or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                or checker.check_reaction("Sonogashira acetylene_aryl OTf", rsmi)
                or checker.check_reaction("Sonogashira alkyne_aryl OTf", rsmi)
                or
                # Buchwald-Hartwig/N-arylation
                checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                )
                or checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                )
                or checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                or
                # Heck reactions
                checker.check_reaction("Heck terminal vinyl", rsmi)
                or checker.check_reaction("Heck_terminal_vinyl", rsmi)
                or
                # Negishi coupling
                checker.check_reaction("Negishi coupling", rsmi)
                or
                # Stille coupling
                checker.check_reaction("Stille reaction_aryl", rsmi)
                or checker.check_reaction("Stille reaction_vinyl", rsmi)
            )

            # Check for other common linear reaction types
            is_linear_reaction = (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                )
                or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                or checker.check_reaction("Reductive amination with ketone", rsmi)
                or checker.check_reaction("Wittig reaction with triphenylphosphorane", rsmi)
                or checker.check_reaction("Wittig with Phosphonium", rsmi)
            )

            if is_coupling_reaction:
                print(f"Coupling reaction detected (considered linear): {rsmi}")
            elif is_linear_reaction:
                print(f"Linear reaction type detected: {rsmi}")
            else:
                # Count significant reactants (excluding small molecules/reagents)
                significant_reactants = 0
                has_building_block = False

                for reactant in reactants:
                    if not reactant:
                        continue

                    mol = Chem.MolFromSmiles(reactant)
                    if not mol:
                        continue

                    # Check for common building blocks
                    is_building_block = (
                        # Boronic compounds
                        checker.check_fg("Boronic acid", reactant)
                        or checker.check_fg("Boronic ester", reactant)
                        or
                        # Halides (common coupling partners)
                        checker.check_fg("Aromatic halide", reactant)
                        or checker.check_fg("Primary halide", reactant)
                        or checker.check_fg("Secondary halide", reactant)
                        or
                        # Organometallics
                        checker.check_fg("Magnesium halide", reactant)
                        or checker.check_fg("Zinc halide", reactant)
                        or
                        # Activated alkenes/alkynes
                        checker.check_fg("Alkyne", reactant)
                        or checker.check_fg("Vinyl", reactant)
                        or
                        # Small functional groups
                        (checker.check_fg("Aldehyde", reactant) and mol.GetNumHeavyAtoms() < 12)
                        or (checker.check_fg("Ketone", reactant) and mol.GetNumHeavyAtoms() < 12)
                        or (
                            checker.check_fg("Carboxylic acid", reactant)
                            and mol.GetNumHeavyAtoms() < 12
                        )
                        or (checker.check_fg("Ester", reactant) and mol.GetNumHeavyAtoms() < 12)
                        or (
                            checker.check_fg("Acyl halide", reactant)
                            and mol.GetNumHeavyAtoms() < 12
                        )
                    )

                    if is_building_block:
                        has_building_block = True
                        print(f"Building block detected: {reactant}")

                    # Consider molecules with more than 10 heavy atoms as significant
                    # unless they're common building blocks with <15 atoms
                    if mol.GetNumHeavyAtoms() > 10 and not (
                        is_building_block and mol.GetNumHeavyAtoms() < 15
                    ):
                        significant_reactants += 1

                # If more than one significant reactant and not involving building blocks,
                # it's likely not a linear synthesis
                if significant_reactants > 1 and not has_building_block:
                    print(f"Non-linear (convergent) step detected in reaction: {rsmi}")
                    print(f"  Number of significant reactants: {significant_reactants}")
                    is_linear = False
                elif significant_reactants > 2:
                    # Even with building blocks, having more than 2 significant reactants
                    # suggests a convergent or multi-component reaction
                    print(f"Multi-component reaction detected (non-linear): {rsmi}")
                    print(f"  Number of significant reactants: {significant_reactants}")
                    is_linear = False

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return is_linear
