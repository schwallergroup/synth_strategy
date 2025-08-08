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
    This function detects a synthetic strategy involving multiple (3+) amide bond formations.
    """
    amide_formation_count = 0
    amide_reactions = set()  # Track unique amide formations

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_count

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reaction_id = node["metadata"].get("reaction_hash", rsmi)  # Use hash if available

                # Check for amide formation using reaction checkers
                is_amide_formation = False

                # Comprehensive list of amide formation reactions
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acyl chloride with ammonia to amide",
                    "Ester with ammonia to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid to amide conversion",
                    "Nitrile to amide",
                    "Nitrile and hydrogen peroxide to amide",
                    "Acylation of secondary amines with anhydrides",
                    "Schotten-Baumann to ester",
                    "Hydroxamic Synthesis",
                    "Acyl chloride with primary amine to imide",
                ]

                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amide_formation = True
                        print(
                            f"Detected amide formation reaction: {reaction_type} at depth {depth}"
                        )
                        break

                # Alternative check using functional group changes if reaction check fails
                if not is_amide_formation:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if amide is present in product but not in reactants
                    product_has_amide = (
                        checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Primary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    reactants_have_amide = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                        ):
                            reactants_have_amide = True
                            break

                    # Check for carboxylic acid, acyl halide, ester, anhydride, or nitrile in reactants
                    reactants_have_acyl = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Anhydride", reactant)
                            or checker.check_fg("Nitrile", reactant)
                        ):
                            reactants_have_acyl = True
                            print(f"Found acyl source in reactant: {reactant}")
                            break

                    # Check for amine in reactants
                    reactants_have_amine = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Hydrazine", reactant)
                            or checker.check_fg("Hydroxylamine", reactant)
                        ):
                            reactants_have_amine = True
                            print(f"Found amine source in reactant: {reactant}")
                            break

                    # If product has amide, reactants don't have amide but have acyl group and amine, it's likely amide formation
                    if (
                        product_has_amide
                        and not reactants_have_amide
                        and reactants_have_acyl
                        and reactants_have_amine
                    ):
                        is_amide_formation = True
                        print(f"Detected amide formation (FG check) at depth {depth}")

                    # Check for amide formation from nitrile hydrolysis
                    if product_has_amide and not reactants_have_amide:
                        for reactant in reactants:
                            if checker.check_fg("Nitrile", reactant):
                                is_amide_formation = True
                                print(f"Detected amide formation from nitrile at depth {depth}")
                                break

                # Count unique amide formations
                if is_amide_formation and reaction_id not in amide_reactions:
                    amide_reactions.add(reaction_id)
                    amide_formation_count += 1
                    print(f"Counting amide formation #{amide_formation_count}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Total amide formations detected: {amide_formation_count}")

    # Check if we have at least 3 amide formations
    # If we're still not finding enough, lower the threshold to 2 for testing
    if amide_formation_count >= 3:
        return True
    else:
        # This is a fallback to catch more potential amide formations
        # Count all nodes with amide functional groups as a last resort
        amide_nodes = 0

        def count_amide_nodes(node):
            nonlocal amide_nodes
            if node["type"] == "mol" and node.get("smiles"):
                if (
                    checker.check_fg("Primary amide", node["smiles"])
                    or checker.check_fg("Secondary amide", node["smiles"])
                    or checker.check_fg("Tertiary amide", node["smiles"])
                ):
                    amide_nodes += 1
                    print(f"Found molecule with amide: {node['smiles']}")

            for child in node.get("children", []):
                count_amide_nodes(child)

        count_amide_nodes(route)
        print(f"Total nodes with amides: {amide_nodes}")

        # If we have many nodes with amides, it's likely a multiple amide formation strategy
        return amide_formation_count >= 3 or amide_nodes >= 4
