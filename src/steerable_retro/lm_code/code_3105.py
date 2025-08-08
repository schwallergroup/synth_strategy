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
    Detects a strategy with early heterocyclic ring formation followed by
    sequential late-stage functionalization (acylation then halogenation)
    """
    # Track if we found the key features
    found_heterocycle = False
    found_acylation = False
    found_halogenation = False
    heterocycle_depth = -1
    acylation_depth = -1
    halogenation_depth = -1

    # Store the final product SMILES
    final_product_smiles = route["smiles"] if route["type"] == "mol" else ""

    # List of heterocyclic rings to check
    heterocycles = [
        "furan",
        "pyrrole",
        "thiophene",
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
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    # Check if the final product contains a heterocycle
    final_product_heterocycle = None
    for heterocycle in heterocycles:
        if checker.check_ring(heterocycle, final_product_smiles):
            final_product_heterocycle = heterocycle
            print(f"Final product contains heterocycle: {heterocycle}")
            found_heterocycle = True
            break

    if not found_heterocycle:
        print("Final product does not contain a heterocycle")
        return False

    # List of acylation reactions to check
    acylation_reactions = [
        "Friedel-Crafts acylation",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of secondary amines",
        "Acylation of primary amines",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    ]

    # List of halogenation reactions to check
    halogenation_reactions = [
        "Aromatic fluorination",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Chlorination",
        "Fluorination",
        "Iodination",
        "Bromination",
    ]

    # List of acyl functional groups to check
    acyl_fgs = [
        "Ketone",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Ester",
        "Carboxylic acid",
        "Acyl halide",
    ]

    # List of halogen functional groups to check
    halogen_fgs = [
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Alkenyl halide",
        "Haloalkyne",
        "Triflate",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle, found_acylation, found_halogenation
        nonlocal heterocycle_depth, acylation_depth, halogenation_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Create RDKit molecules
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                # Skip if any molecule failed to parse
                if product_mol is None or None in reactant_mols:
                    print(f"Failed to parse molecules at depth {depth}")
                    return

                # Check for heterocyclic ring formation
                if not found_heterocycle or heterocycle_depth == -1:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product_smiles):
                            # Verify this heterocycle wasn't already in the reactants
                            if not any(
                                checker.check_ring(heterocycle, r) for r in reactants_smiles
                            ):
                                found_heterocycle = True
                                heterocycle_depth = depth
                                print(
                                    f"Found heterocyclic ring formation ({heterocycle}) at depth {depth}"
                                )
                                break

                # Check for acylation reactions
                if not found_acylation:
                    for acylation_reaction in acylation_reactions:
                        if checker.check_reaction(acylation_reaction, rsmi):
                            # Verify the acylation occurs on a molecule containing the heterocycle
                            if checker.check_ring(final_product_heterocycle, product_smiles):
                                found_acylation = True
                                acylation_depth = depth
                                print(
                                    f"Found acylation reaction ({acylation_reaction}) at depth {depth}"
                                )
                                break

                # If no specific acylation reaction was found, check for acyl group addition
                if not found_acylation:
                    # Check for acyl functional groups in product that weren't in reactants
                    for acyl_fg in acyl_fgs:
                        if checker.check_fg(acyl_fg, product_smiles):
                            if not any(checker.check_fg(acyl_fg, r) for r in reactants_smiles):
                                # Verify the molecule with the new acyl group contains the heterocycle
                                if checker.check_ring(final_product_heterocycle, product_smiles):
                                    found_acylation = True
                                    acylation_depth = depth
                                    print(f"Found acylation (new {acyl_fg}) at depth {depth}")
                                    break

                # Check for halogenation reactions
                if not found_halogenation:
                    for halogenation_reaction in halogenation_reactions:
                        if checker.check_reaction(halogenation_reaction, rsmi):
                            # Verify the halogenation occurs on a molecule containing the heterocycle
                            if checker.check_ring(final_product_heterocycle, product_smiles):
                                found_halogenation = True
                                halogenation_depth = depth
                                print(
                                    f"Found halogenation reaction ({halogenation_reaction}) at depth {depth}"
                                )
                                break

                # If no specific halogenation reaction was found, check for halogen addition
                if not found_halogenation:
                    # Check for halogen functional groups in product that weren't in reactants
                    for halogen_fg in halogen_fgs:
                        if checker.check_fg(halogen_fg, product_smiles):
                            if not any(checker.check_fg(halogen_fg, r) for r in reactants_smiles):
                                # Verify the molecule with the new halogen group contains the heterocycle
                                if checker.check_ring(final_product_heterocycle, product_smiles):
                                    found_halogenation = True
                                    halogenation_depth = depth
                                    print(f"Found halogenation (new {halogen_fg}) at depth {depth}")
                                    break

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Summary - Heterocycle: {found_heterocycle} at depth {heterocycle_depth}, "
        f"Acylation: {found_acylation} at depth {acylation_depth}, "
        f"Halogenation: {found_halogenation} at depth {halogenation_depth}"
    )

    # Check if we found the pattern
    if found_heterocycle and found_acylation and found_halogenation:
        # Case 1: We found explicit heterocycle formation in the route
        if heterocycle_depth != -1:
            # Check proper sequence: heterocycle formation (early) -> acylation -> halogenation (late)
            if heterocycle_depth > acylation_depth > halogenation_depth:
                print(
                    "Found strategy: early heterocycle formation with late-stage functionalization (explicit formation)"
                )
                return True

        # Case 2: Heterocycle was present from the beginning, check if functionalization follows proper sequence
        elif acylation_depth > halogenation_depth:
            print(
                "Found strategy: heterocycle present from beginning with late-stage functionalization"
            )
            return True

        # Case 3: Alternative order where halogenation is definitely the last step
        elif halogenation_depth < acylation_depth:
            print("Found strategy: heterocycle with late-stage halogenation as final step")
            return True

    return False
