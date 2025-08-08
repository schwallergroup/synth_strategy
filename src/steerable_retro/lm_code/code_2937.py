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
    Detects if the synthetic route involves a sequential SNAr strategy with
    multiple nucleophilic aromatic substitutions.
    """
    snar_reactions = []

    def calculate_depth(current_node, root_node, current_depth=0):
        """Helper function to calculate depth if not provided in metadata"""
        if current_node == root_node:
            return 0

        # Traverse the tree to find the current node's depth
        def find_node_depth(node, target, depth=0):
            if node == target:
                return depth

            for child in node.get("children", []):
                result = find_node_depth(child, target, depth + 1)
                if result is not None:
                    return result
            return None

        return find_node_depth(root_node, current_node)

    def dfs_traverse(node, current_depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Get depth from metadata or calculate it
            if "depth" in node["metadata"]:
                depth = int(node["metadata"]["depth"])
            else:
                depth = current_depth

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for nucleophilic aromatic substitution reactions using predefined reactions
            is_snar = (
                checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("N-arylation", rsmi)
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution thiol", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution aryl alcohol", rsmi)
                or checker.check_reaction("Ullmann condensation", rsmi)
            )

            if is_snar:
                print(f"Found SNAr reaction through reaction check at depth {depth}")
                snar_reactions.append(depth)
            else:
                # If not found in predefined reactions, check for the pattern manually
                print("Checking for SNAr pattern manually...")

                # Check if any reactant has an aromatic ring with a leaving group
                aromatic_with_leaving_group = False
                nucleophile_present = False
                has_ewg = False

                # Define common aromatic rings
                aromatic_rings = [
                    "benzene",
                    "pyridine",
                    "pyrimidine",
                    "pyrazine",
                    "pyridazine",
                    "furan",
                    "thiophene",
                    "pyrrole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "triazole",
                    "tetrazole",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                ]

                for reactant in reactants:
                    try:
                        # Check for aromatic rings
                        has_aromatic_ring = any(
                            checker.check_ring(ring, reactant) for ring in aromatic_rings
                        )

                        if has_aromatic_ring:
                            print(f"Found aromatic ring in reactant: {reactant}")

                            # Check for leaving groups attached to aromatic rings
                            if (
                                checker.check_fg("Aromatic halide", reactant)
                                or checker.check_fg("Triflate", reactant)
                                or checker.check_fg("Tosylate", reactant)
                                or checker.check_fg("Mesylate", reactant)
                                or checker.check_fg("Thiocyanate", reactant)
                                or "SCN" in reactant
                                or "SC#N" in reactant
                            ):
                                print(f"Found leaving group in reactant: {reactant}")
                                aromatic_with_leaving_group = True

                            # Check for electron-withdrawing groups that activate SNAr
                            if (
                                checker.check_fg("Nitro group", reactant)
                                or checker.check_fg("Nitrile", reactant)
                                or checker.check_fg("Ester", reactant)
                                or checker.check_fg("Ketone", reactant)
                                or checker.check_fg("Trifluoro group", reactant)
                                or checker.check_fg("Sulfone", reactant)
                                or checker.check_fg("Sulfonamide", reactant)
                            ):
                                print(f"Found electron-withdrawing group in reactant: {reactant}")
                                has_ewg = True

                        # Check for nucleophiles
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Phenol", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Aromatic thiol", reactant)
                            or checker.check_fg("Aliphatic thiol", reactant)
                            or checker.check_fg("Aniline", reactant)
                            or "SH" in reactant
                        ):
                            print(f"Found nucleophile in reactant: {reactant}")
                            nucleophile_present = True
                    except Exception as e:
                        print(f"Error analyzing reactant: {e}")
                        continue

                # Verify the product has an aromatic ring with the nucleophile attached
                product_has_substitution = False
                try:
                    # Check if product has expected substitution patterns
                    has_product_aromatic_ring = any(
                        checker.check_ring(ring, product) for ring in aromatic_rings
                    )

                    if has_product_aromatic_ring:
                        print(f"Found aromatic ring in product: {product}")
                        if (
                            checker.check_fg("Aniline", product)
                            or checker.check_fg("Phenol", product)
                            or checker.check_fg("Ether", product)
                            or checker.check_fg("Monosulfide", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        ):
                            print(f"Found nucleophilic substitution in product: {product}")
                            product_has_substitution = True
                except Exception as e:
                    print(f"Error analyzing product: {e}")

                # For heterocycles, we might not need an EWG
                is_heterocycle = any(
                    checker.check_ring(ring, reactant)
                    for ring in [
                        "pyridine",
                        "pyrimidine",
                        "pyrazine",
                        "pyridazine",
                        "triazole",
                        "tetrazole",
                        "oxazole",
                        "thiazole",
                    ]
                    for reactant in reactants
                )

                # Special case for thiocyanate displacement
                thiocyanate_displacement = False
                for reactant in reactants:
                    if "SC#N" in reactant or "SCN" in reactant:
                        print(f"Found thiocyanate in reactant: {reactant}")
                        thiocyanate_displacement = True

                # SNAr conditions:
                # 1. Has aromatic ring with leaving group
                # 2. Has nucleophile
                # 3. Product shows substitution
                # 4. Either has EWG or is a heterocycle (which can be activated without EWG)
                manual_snar = (
                    (aromatic_with_leaving_group or thiocyanate_displacement)
                    and nucleophile_present
                    and (product_has_substitution or is_heterocycle)
                    and (has_ewg or is_heterocycle)
                )

                if manual_snar:
                    print(f"Found SNAr reaction through manual check at depth {depth}")
                    snar_reactions.append(depth)

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    print("Starting traversal of synthesis route")
    dfs_traverse(route)

    # Check if there are multiple SNAr reactions at different depths
    if len(snar_reactions) >= 2:
        # Sort reactions by depth to check if they form a sequence
        snar_reactions.sort()
        print(f"Found SNAr reactions at depths: {snar_reactions}")

        # Check if the reactions are at different depths (sequential)
        if len(set(snar_reactions)) >= 2:
            print("Confirmed sequential SNAr strategy")
            return True
        else:
            print("Found multiple SNAr reactions but they are at the same depth")
    else:
        print(
            f"Found only {len(snar_reactions)} SNAr reactions, need at least 2 for a sequential strategy"
        )

    return False
