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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    Linear synthesis is characterized by each reaction having only one non-reagent reactant.
    """
    is_linear = True
    reaction_count = 0
    convergent_step_count = 0

    # Common reagents and reagent-like functional groups
    reagent_fg_types = [
        "Triflate",
        "Mesylate",
        "Tosylate",
        "Magnesium halide",
        "Zinc halide",
        "Tin",
        "Alkyl lithium",
        "Aryl lithium",
        "Silane",
        "Silyl protective group",
        "TMS ether protective group",
        "Boronic acid",
        "Boronic ester",
        "Boc",
        "Acetal/Ketal",
        "Carbon dioxide",
        "Carbon monoxide",
        "Methanol",
        "Oxime",
    ]

    # Small ring structures that indicate building blocks rather than reagents
    building_block_rings = [
        "benzene",
        "pyridine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
    ]

    # Reactions that are typically convergent
    convergent_reaction_types = [
        "Suzuki coupling",
        "Negishi coupling",
        "Stille reaction",
        "Heck reaction",
        "Sonogashira",
        "Buchwald-Hartwig",
        "Ullmann condensation",
        "Ugi reaction",
        "A3 coupling",
        "Petasis reaction",
        "Chan-Lam",
        "Kumada cross-coupling",
        "Hiyama-Denmark Coupling",
        "Mitsunobu",
        "Aldol condensation",
        "Wittig reaction",
        "Diels-Alder",
        "Michael addition",
        "Reductive amination",
    ]

    def is_likely_reagent(smiles):
        """Determine if a molecule is likely a reagent rather than a key reactant."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False

            # Very small molecules are often reagents, unless they contain important ring structures
            if mol.GetNumAtoms() <= 6:
                # Check if it contains any building block rings
                if any(checker.check_ring(ring, smiles) for ring in building_block_rings):
                    return False
                return True

            # Check for reagent-like functional groups
            for fg in reagent_fg_types:
                if checker.check_fg(fg, smiles):
                    return True

            # Common reagent patterns
            if "[Mg]" in smiles or "[Li]" in smiles or "[Zn]" in smiles:  # Organometallics
                return True
            if "Si(" in smiles:  # Silyl compounds
                return True

            return False
        except Exception as e:
            print(f"Error in is_likely_reagent: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count, convergent_step_count

        if node["type"] == "reaction":
            reaction_count += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a known convergent reaction type
                is_convergent_rxn_type = False
                for rxn_type in convergent_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found convergent reaction type '{rxn_type}' at depth {depth}")
                        is_convergent_rxn_type = True
                        convergent_step_count += 1
                        break

                # If not already identified as convergent by reaction type, check reactants
                if not is_convergent_rxn_type:
                    # Count non-reagent reactants
                    non_reagent_count = 0
                    for reactant in reactants:
                        if not is_likely_reagent(reactant):
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.GetNumAtoms() > 3:  # More conservative threshold
                                non_reagent_count += 1

                    # If more than one significant reactant, it's potentially a convergent step
                    if non_reagent_count > 1:
                        print(
                            f"Found convergent step at depth {depth} with {non_reagent_count} significant reactants"
                        )
                        convergent_step_count += 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Allow for a small percentage of convergent steps in an otherwise linear synthesis
    # Only consider it a linear synthesis if we have at least 2 reactions
    if reaction_count >= 2:
        # If more than 25% of steps are convergent, consider it non-linear
        convergent_percentage = convergent_step_count / reaction_count
        print(
            f"Convergent percentage: {convergent_percentage:.2f} ({convergent_step_count}/{reaction_count})"
        )
        if convergent_percentage > 0.25:
            is_linear = False
        print(
            f"Synthesis has {reaction_count} reactions with {convergent_step_count} convergent steps"
        )
        print(f"Final decision: {'Linear' if is_linear else 'Convergent'} synthesis")
        return is_linear
    else:
        print(
            f"Synthesis has only {reaction_count} reactions, which is insufficient to determine strategy"
        )
        return False  # Not enough reactions to determine strategy
