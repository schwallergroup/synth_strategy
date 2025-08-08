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
    Detects a linear synthesis strategy where a heterocyclic core is elaborated through
    sequential functionalization steps rather than convergent assembly.
    """
    # Count the number of reactions
    reaction_count = 0
    # Count reactions with multiple reactants (potential convergent steps)
    convergent_steps = 0
    # Track heterocycle presence and when it first appears
    heterocycle_nodes = []
    # Track elaboration steps (reactions that modify the heterocycle)
    elaboration_steps = 0
    # Track max depth of the synthesis
    max_depth = 0
    # Track sequential modifications of the same heterocycle
    sequential_modifications = 0
    # Previous heterocycle type
    prev_heterocycle_type = None

    # List of heterocycles to check
    heterocycle_types = [
        "furan",
        "pyran",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, convergent_steps, heterocycle_nodes, elaboration_steps, max_depth
        nonlocal sequential_modifications, prev_heterocycle_type

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            # Check for heterocycles in molecule nodes
            mol_smiles = node["smiles"]
            for heterocycle in heterocycle_types:
                if checker.check_ring(heterocycle, mol_smiles):
                    heterocycle_nodes.append((mol_smiles, depth, heterocycle))
                    print(f"Found heterocycle {heterocycle} at depth {depth}: {mol_smiles}")

                    # Check if this is the same heterocycle type as previous
                    if prev_heterocycle_type == heterocycle:
                        sequential_modifications += 1
                        print(f"Sequential modification of {heterocycle} detected")

                    prev_heterocycle_type = heterocycle
                    break

        elif node["type"] == "reaction":
            reaction_count += 1

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a convergent step (multiple reactants)
                if len(reactants) > 1:
                    convergent_steps += 1
                    print(f"Found convergent step with {len(reactants)} reactants: {rsmi}")

                # Check if this reaction modifies a heterocycle
                # First, check if product contains a heterocycle
                product_has_heterocycle = False
                product_heterocycle_type = None
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, product):
                        product_has_heterocycle = True
                        product_heterocycle_type = heterocycle
                        break

                # Then check if any reactant contains the same heterocycle
                if product_has_heterocycle:
                    for reactant in reactants:
                        if checker.check_ring(product_heterocycle_type, reactant):
                            # This reaction maintains the heterocycle, check if it's an elaboration
                            # Look for common functionalization reactions
                            elaboration_reaction_types = [
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                "N-alkylation of primary amines with alkyl halides",
                                "N-alkylation of secondary amines with alkyl halides",
                                "Suzuki coupling with boronic acids",
                                "Suzuki coupling with boronic esters",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                "Friedel-Crafts acylation",
                                "Friedel-Crafts alkylation",
                                "Aromatic nitration with HNO3",
                                "Aromatic fluorination",
                                "Aromatic chlorination",
                                "Aromatic bromination",
                                "Aromatic iodination",
                            ]

                            for rxn_type in elaboration_reaction_types:
                                if checker.check_reaction(rxn_type, rsmi):
                                    elaboration_steps += 1
                                    print(f"Found elaboration step: {rxn_type} at depth {depth}")
                                    break

                            # Also check for general functional group changes on the heterocycle
                            p_mol = Chem.MolFromSmiles(product)
                            r_mol = Chem.MolFromSmiles(reactant)
                            if p_mol and r_mol:
                                # Check if functional groups have changed
                                fg_types = [
                                    "Primary amine",
                                    "Secondary amine",
                                    "Tertiary amine",
                                    "Primary alcohol",
                                    "Secondary alcohol",
                                    "Tertiary alcohol",
                                    "Carboxylic acid",
                                    "Ester",
                                    "Amide",
                                    "Nitrile",
                                    "Primary halide",
                                    "Secondary halide",
                                    "Tertiary halide",
                                    "Aromatic halide",
                                ]

                                for fg in fg_types:
                                    # Check if the functional group status changed
                                    if checker.check_fg(fg, reactant) != checker.check_fg(
                                        fg, product
                                    ):
                                        # Verify the change is related to the heterocycle
                                        # This is a simplification - in a full implementation we would
                                        # use atom mapping to verify the change is on/near the heterocycle
                                        elaboration_steps += 1
                                        print(
                                            f"Found functional group change: {fg} at depth {depth}"
                                        )
                                        break

                            break

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Calculate the ratio of convergent steps to total reactions
    convergent_ratio = convergent_steps / reaction_count if reaction_count > 0 else 0

    # Calculate percentage of elaboration steps relative to total reactions
    elaboration_percentage = elaboration_steps / reaction_count if reaction_count > 0 else 0

    # Find the earliest heterocycle (lowest depth)
    earliest_heterocycle_depth = 0
    if heterocycle_nodes:
        earliest_heterocycle_depth = min([depth for _, depth, _ in heterocycle_nodes])
        print(
            f"Earliest heterocycle found at depth {earliest_heterocycle_depth} (max depth: {max_depth})"
        )

    # Determine if heterocycle appears early in the synthesis (in the first half)
    heterocycle_appears_early = (
        earliest_heterocycle_depth < max_depth / 2 if max_depth > 0 else False
    )

    print(f"Reaction count: {reaction_count}")
    print(f"Convergent steps: {convergent_steps}")
    print(f"Elaboration steps: {elaboration_steps}")
    print(f"Sequential modifications: {sequential_modifications}")
    print(f"Heterocycle appears early: {heterocycle_appears_early}")
    print(f"Convergent ratio: {convergent_ratio}")
    print(f"Elaboration percentage: {elaboration_percentage}")

    # If the route has multiple steps, contains heterocycles early, has elaboration steps, and is primarily linear
    # or has a significant proportion of elaboration steps
    if (
        reaction_count >= 3
        and heterocycle_nodes
        and heterocycle_appears_early
        and elaboration_steps >= 1
        and (convergent_ratio <= 0.7 or elaboration_percentage >= 0.5)
    ):
        print(
            f"Found linear heterocycle elaboration strategy with {convergent_steps}/{reaction_count} convergent steps"
        )
        print(f"Heterocycle elaboration steps: {elaboration_steps}")
        return True

    return False
