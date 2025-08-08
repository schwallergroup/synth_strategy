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
    Detects a synthesis strategy involving sequential nucleophilic substitutions at a benzylic position.
    """
    # List of nucleophilic substitution reaction types to check
    nucleophilic_rxn_types = [
        "Williamson Ether Synthesis",
        "S-alkylation of thiols",
        "S-alkylation of thiols (ethyl)",
        "S-alkylation of thiols with alcohols",
        "S-alkylation of thiols with alcohols (ethyl)",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Alcohol to chloride_sulfonyl chloride",
        "Alcohol to chloride_SOCl2",
        "Alcohol to chloride_HCl",
        "Alcohol to chloride_Salt",
        "Alcohol to chloride_Other",
        "Finkelstein reaction",
        "thioether_nucl_sub",
        "Appel reaction",
        "Mitsunobu aryl ether",
        "Mitsunobu esterification",
        "Mitsunobu aryl ether (intramolecular)",
        "Mitsunobu_phenole",
        "Mitsunobu_imide",
        "Mitsunobu_sulfonamide",
        "Mitsunobu_tetrazole_1",
        "Mitsunobu_tetrazole_2",
        "Mitsunobu_tetrazole_3",
        "Mitsunobu_tetrazole_4",
        "Alcohol to azide",
        "Amine to azide",
    ]

    # List of leaving groups to check
    leaving_groups = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Triflate",
        "Mesylate",
        "Tosylate",
        "Primary alcohol",  # Can be a leaving group in some reactions (e.g., Mitsunobu)
        "Secondary alcohol",
        "Tertiary alcohol",
    ]

    # List of nucleophilic groups
    nucleophiles = [
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Aliphatic thiol",
        "Phenol",
        "Azide",
        "Carboxylic acid",
    ]

    # Helper function to detect nucleophilic substitution
    def is_nucleophilic_substitution(rsmi):
        # Method 1: Check if it matches known nucleophilic substitution reaction types
        for rxn_type in nucleophilic_rxn_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected nucleophilic substitution: {rxn_type}")
                return True

        # Method 2: Check for benzylic substitution pattern
        reactants_smiles = rsmi.split(">")[0].split(".")
        product_smiles = rsmi.split(">")[-1]

        # Check if either product or any reactant has a benzene ring
        product_has_benzene = checker.check_ring("benzene", product_smiles)

        if product_has_benzene:
            # Check if any reactant has a leaving group and another has a nucleophile
            has_leaving_group = False
            has_nucleophile = False

            for reactant in reactants_smiles:
                if not reactant:
                    continue

                # Check for leaving groups
                for leaving_group in leaving_groups:
                    if checker.check_fg(leaving_group, reactant):
                        has_leaving_group = True
                        break

                # Check for nucleophiles
                for nucleophile in nucleophiles:
                    if checker.check_fg(nucleophile, reactant):
                        has_nucleophile = True
                        break

            # If we have both a leaving group and a nucleophile in the reactants,
            # and the product has a benzene ring, it's likely a benzylic substitution
            if has_leaving_group and has_nucleophile:
                print(
                    f"Detected benzylic nucleophilic substitution based on functional group patterns"
                )
                return True

        # Check reactants for benzene rings
        for reactant in reactants_smiles:
            if not reactant:
                continue

            if checker.check_ring("benzene", reactant):
                # Check if this reactant has a leaving group
                has_leaving_group = any(checker.check_fg(lg, reactant) for lg in leaving_groups)

                # Check if any other reactant has a nucleophile
                has_nucleophile = False
                for other_reactant in reactants_smiles:
                    if other_reactant != reactant:
                        has_nucleophile = any(
                            checker.check_fg(nuc, other_reactant) for nuc in nucleophiles
                        )
                        if has_nucleophile:
                            break

                if has_leaving_group and has_nucleophile:
                    print(f"Detected benzylic nucleophilic substitution (reactant has benzene)")
                    return True

        return False

    # Track reaction nodes with nucleophilic substitutions
    substitution_reactions = []

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        # Create a copy of the current path and add this node
        current_path = path.copy()
        current_path.append(node)

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a nucleophilic substitution reaction
            if is_nucleophilic_substitution(rsmi):
                substitution_reactions.append((depth, node, current_path))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Sort substitution reactions by depth (to find sequential reactions)
    substitution_reactions.sort(key=lambda x: x[0])

    print(f"Found {len(substitution_reactions)} nucleophilic substitution reactions")

    # Check if we have at least 2 nucleophilic substitutions
    if len(substitution_reactions) >= 2:
        # Check if any two substitutions are in the same path
        for i in range(len(substitution_reactions) - 1):
            for j in range(i + 1, len(substitution_reactions)):
                depth_i, node_i, path_i = substitution_reactions[i]
                depth_j, node_j, path_j = substitution_reactions[j]

                # Check if the paths share common molecules (indicating sequential reactions)
                # We need to extract molecule nodes from the paths
                mol_nodes_i = [n for n in path_i if n["type"] == "mol"]
                mol_nodes_j = [n for n in path_j if n["type"] == "mol"]

                # Check for common molecules
                common_mols = set(n["smiles"] for n in mol_nodes_i) & set(
                    n["smiles"] for n in mol_nodes_j
                )

                if common_mols:
                    print(
                        f"Found sequential nucleophilic substitutions with common molecules: {common_mols}"
                    )
                    return True

    # If no sequential nucleophilic substitutions found, return False
    return False
