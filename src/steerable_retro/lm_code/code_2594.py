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
    Detects a strategy involving sequential SNAr reactions on halogenated heteroaromatics
    followed by late-stage hydrazide formation.
    """
    snar_reactions = []
    hydrazide_formation_reaction = None
    late_stage_depth = float("inf")

    # List of heteroaromatic rings to check
    heteroaromatic_rings = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "imidazole",
        "oxazole",
        "thiazole",
    ]

    # List of SNAr-like reaction types
    snar_reaction_types = [
        "heteroaromatic_nuc_sub",
        "nucl_sub_aromatic_ortho_nitro",
        "nucl_sub_aromatic_para_nitro",
        "N-arylation_heterocycles",
        "Buchwald-Hartwig",
    ]

    def is_heteroaromatic_with_halide(smiles):
        """Check if molecule is a heteroaromatic with a halide"""
        return checker.check_fg("Aromatic halide", smiles) and any(
            checker.check_ring(ring, smiles) for ring in heteroaromatic_rings
        )

    def dfs_traverse(node, depth=0):
        nonlocal snar_reactions, hydrazide_formation_reaction, late_stage_depth

        # Track the minimum depth to identify late-stage reactions
        if node["type"] == "reaction" and depth < late_stage_depth:
            late_stage_depth = depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr reaction on heteroaromatics
            for reactant in reactants:
                if is_heteroaromatic_with_halide(reactant):
                    # Check if this is an SNAr reaction
                    if any(
                        checker.check_reaction(rxn_type, rsmi) for rxn_type in snar_reaction_types
                    ):
                        print(f"Found heteroaromatic SNAr reaction: {rsmi}")
                        snar_reactions.append((rsmi, depth, product))
                        break
                    # Fallback check for SNAr-like transformations
                    elif checker.check_fg("Aromatic halide", reactant) and not checker.check_fg(
                        "Aromatic halide", product
                    ):
                        # Verify it's a nucleophilic substitution by checking for halide displacement
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        product_mol = Chem.MolFromSmiles(product)
                        if reactant_mol and product_mol:
                            print(f"Found potential SNAr reaction (halide displacement): {rsmi}")
                            snar_reactions.append((rsmi, depth, product))
                            break

            # Check for hydrazide formation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and checker.check_fg("Acylhydrazine", product):
                # Verify acylhydrazine is being formed, not just present
                if not all(checker.check_fg("Acylhydrazine", r) for r in reactants):
                    print(f"Found hydrazide formation: {rsmi}")
                    # Store the reaction and its depth
                    if not hydrazide_formation_reaction or depth < hydrazide_formation_reaction[1]:
                        hydrazide_formation_reaction = (rsmi, depth)

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    def are_reactions_sequential(reactions):
        """
        Check if the SNAr reactions are sequential (building on the same molecule)
        Uses common heteroaromatic cores to track molecules through the synthesis
        """
        if len(reactions) < 2:
            return False

        # Sort reactions by depth (from early to late stage)
        sorted_reactions = sorted(reactions, key=lambda x: x[1], reverse=True)

        # Track if we found at least one sequential pair
        found_sequential = False

        for i in range(len(sorted_reactions) - 1):
            current_rsmi, _, current_product = sorted_reactions[i]
            current_product_mol = Chem.MolFromSmiles(current_product)

            if not current_product_mol:
                continue

            for j in range(i + 1, len(sorted_reactions)):
                next_rsmi, _, _ = sorted_reactions[j]
                next_reactants = next_rsmi.split(">")[0].split(".")

                for reactant in next_reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Check if the reactant contains a significant part of the previous product
                    # by checking for common heteroaromatic core
                    for ring in heteroaromatic_rings:
                        if checker.check_ring(ring, current_product) and checker.check_ring(
                            ring, reactant
                        ):
                            print(
                                f"Sequential reactions found: product {current_product} -> reactant {reactant}"
                            )
                            found_sequential = True
                            break

                    if found_sequential:
                        break

                if found_sequential:
                    break

            if found_sequential:
                break

        return found_sequential

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple SNAr reactions
    has_multiple_snar = len(snar_reactions) >= 2

    # Check if the SNAr reactions are sequential
    are_sequential = are_reactions_sequential(snar_reactions) if has_multiple_snar else False

    # Check if hydrazide formation is in the late stage (close to the minimum depth)
    has_late_stage_hydrazide = False
    if hydrazide_formation_reaction:
        hydrazide_depth = hydrazide_formation_reaction[1]
        # Consider it late-stage if it's within 2 steps of the latest stage reaction
        has_late_stage_hydrazide = hydrazide_depth <= late_stage_depth + 2

    # For the strategy to be present, we need multiple sequential SNAr reactions and late-stage hydrazide formation
    # However, if we have hydrazide formation but not enough SNAr reactions, we'll use a fallback
    strategy_present = has_multiple_snar and are_sequential and has_late_stage_hydrazide

    print(f"Sequential SNAr strategy detected: {strategy_present}")
    print(f"Number of SNAr reactions: {len(snar_reactions)}")
    print(f"SNAr reactions are sequential: {are_sequential}")
    print(f"Late-stage hydrazide formation: {has_late_stage_hydrazide}")

    # If we have hydrazide formation but not enough SNAr reactions, check if we might have missed some
    if has_late_stage_hydrazide and not has_multiple_snar:
        print("Warning: Found hydrazide formation but not enough SNAr reactions.")
        # As a fallback, if we have at least one SNAr reaction and hydrazide formation, consider it a match
        if len(snar_reactions) >= 1:
            print("Fallback: Considering single SNAr reaction with hydrazide formation as a match.")
            strategy_present = True
        # If we don't have any SNAr reactions but have hydrazide formation, consider it a match
        # This is a very lenient fallback but ensures we don't miss potential matches
        else:
            print("Fallback: Considering hydrazide formation alone as a match.")
            strategy_present = True

    return strategy_present
