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
    Detects if the synthesis route involves late-stage coupling of two or more
    complex fragments in the final steps.
    """
    has_late_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_coupling

        if node["type"] == "reaction" and depth <= 2:  # Final, penultimate, or antepenultimate step
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Only consider reactions with multiple valid reactants
                if len(reactants_smiles) >= 2 and all(len(r.strip()) > 0 for r in reactants_smiles):
                    # Check if this is a known coupling reaction
                    is_coupling_reaction = (
                        checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                        or checker.check_reaction("Negishi coupling", rsmi)
                        or checker.check_reaction("Stille reaction_aryl", rsmi)
                        or checker.check_reaction("Stille reaction_vinyl", rsmi)
                        or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        )
                        or checker.check_reaction("Heck terminal vinyl", rsmi)
                        or checker.check_reaction("Ullmann condensation", rsmi)
                        or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                        or checker.check_reaction("Kumada cross-coupling", rsmi)
                    )

                    if is_coupling_reaction:
                        reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                        product = Chem.MolFromSmiles(product_smiles)

                        if product and all(r for r in reactants):
                            # Check complexity of reactants (using number of heavy atoms and excluding common reagents)
                            complex_fragments = 0
                            for r in reactants:
                                r_smiles = Chem.MolToSmiles(r)
                                if (
                                    r.GetNumHeavyAtoms() >= 8
                                    and not checker.check_fg("Magnesium halide", r_smiles)
                                    and not checker.check_fg("Zinc halide", r_smiles)
                                    and not checker.check_fg("Alkyl lithium", r_smiles)
                                    and not checker.check_fg("Boronic acid", r_smiles)
                                    and not checker.check_fg("Boronic ester", r_smiles)
                                ):
                                    complex_fragments += 1

                            if complex_fragments >= 2:
                                print(
                                    f"Late-stage coupling of complex fragments detected at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                                has_late_coupling = True
                    else:
                        # For reactions not in our predefined list, check for C-C bond formation
                        # between two complex fragments
                        reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                        product = Chem.MolFromSmiles(product_smiles)

                        if product and all(r for r in reactants):
                            # Check complexity of reactants
                            complex_reactants = []
                            for r in reactants:
                                r_smiles = Chem.MolToSmiles(r)
                                if (
                                    r.GetNumHeavyAtoms() >= 8
                                    and not checker.check_fg("Magnesium halide", r_smiles)
                                    and not checker.check_fg("Zinc halide", r_smiles)
                                    and not checker.check_fg("Alkyl lithium", r_smiles)
                                ):
                                    complex_reactants.append(r)

                            if len(complex_reactants) >= 2:
                                # This is a potential coupling of complex fragments
                                # Additional check could be performed here to verify bond formation
                                print(
                                    f"Potential late-stage coupling of complex fragments detected at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                                has_late_coupling = True
            except Exception as e:
                print(f"Error processing reaction SMILES at depth {depth}: {str(e)}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_late_coupling
