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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects a convergent synthesis strategy where a heterocyclic fragment
    is combined with a complex ring system in a late-stage reaction.
    """
    # Track if we found the convergent synthesis pattern
    found_convergent_synthesis = False

    # List of heterocycles to check
    heterocycles = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
    ]

    # List of coupling reactions to check
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Heck terminal vinyl",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Sonogashira alkyne_aryl halide",
        "Ullmann condensation",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Goldberg coupling",
        "Stille reaction_vinyl",
        "Sonogashira acetylene_aryl halide",
    ]

    # Track the depth during traversal
    current_depth = 0
    max_depth_to_check = 2  # Check reactions at depth 0, 1, and 2

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_synthesis

        # Check if this is a reaction node with required metadata
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}", flush=True)

            # Only process late-stage reactions (depth 0, 1, or 2)
            if depth <= max_depth_to_check:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")

                # Check if we have at least 2 reactants (convergent)
                if len(reactants) >= 2:
                    print(f"Found convergent reaction with {len(reactants)} reactants", flush=True)

                    # Check for coupling reaction
                    is_coupling = False
                    for reaction_type in coupling_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Identified as {reaction_type}", flush=True)
                            is_coupling = True
                            break

                    # If not a known coupling reaction, check if it's combining fragments
                    if not is_coupling:
                        try:
                            # Check if product has more rings than any individual reactant
                            product_mol = Chem.MolFromSmiles(product_part)
                            if product_mol:
                                product_ring_count = product_mol.GetRingInfo().NumRings()
                                max_reactant_rings = 0
                                for reactant in reactants:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol:
                                        reactant_ring_count = reactant_mol.GetRingInfo().NumRings()
                                        max_reactant_rings = max(
                                            max_reactant_rings, reactant_ring_count
                                        )

                                if product_ring_count > max_reactant_rings:
                                    print(
                                        f"Product has more rings ({product_ring_count}) than any reactant ({max_reactant_rings})",
                                        flush=True,
                                    )
                                    is_coupling = True

                                # Also check for C-C, C-N, or C-O bond formation
                                if not is_coupling:
                                    # Check for common functional groups involved in coupling
                                    for reactant in reactants:
                                        if (
                                            checker.check_fg("Boronic acid", reactant)
                                            or checker.check_fg("Boronic ester", reactant)
                                            or checker.check_fg("Aromatic halide", reactant)
                                            or checker.check_fg("Magnesium halide", reactant)
                                            or checker.check_fg("Zinc halide", reactant)
                                        ):
                                            print(
                                                f"Found coupling-related functional group in {reactant}",
                                                flush=True,
                                            )
                                            is_coupling = True
                                            break
                        except Exception as e:
                            print(f"Error analyzing rings: {e}", flush=True)

                    if is_coupling:
                        # Check for heterocycle in one reactant
                        has_heterocycle = False
                        heterocycle_reactant = None
                        heterocycle_found = None

                        for reactant in reactants:
                            for heterocycle in heterocycles:
                                if checker.check_ring(heterocycle, reactant):
                                    print(
                                        f"Found heterocycle: {heterocycle} in {reactant}",
                                        flush=True,
                                    )
                                    has_heterocycle = True
                                    heterocycle_reactant = reactant
                                    heterocycle_found = heterocycle
                                    break
                            if has_heterocycle:
                                break

                        # Check for complex ring system in another reactant
                        has_complex_ring = False
                        complex_ring_reactant = None

                        for reactant in reactants:
                            # Skip the heterocycle reactant
                            if reactant == heterocycle_reactant:
                                continue

                            try:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    # Consider a complex ring system as having at least 2 rings
                                    ring_count = reactant_mol.GetRingInfo().NumRings()
                                    if ring_count >= 2:
                                        print(
                                            f"Found complex ring system with {ring_count} rings in {reactant}",
                                            flush=True,
                                        )
                                        has_complex_ring = True
                                        complex_ring_reactant = reactant
                                        break
                            except Exception as e:
                                print(f"Error analyzing complex rings: {e}", flush=True)

                        if has_heterocycle and has_complex_ring:
                            print(
                                f"FOUND CONVERGENT SYNTHESIS WITH HETEROCYCLE: {heterocycle_found}",
                                flush=True,
                            )
                            print(f"Heterocycle reactant: {heterocycle_reactant}", flush=True)
                            print(f"Complex ring reactant: {complex_ring_reactant}", flush=True)
                            found_convergent_synthesis = True

        # Continue traversing with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_convergent_synthesis
