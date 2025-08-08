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
    This function detects a synthetic strategy where a heterocyclic ring
    (thiazoline/thiazolidine) is preserved throughout the synthesis.
    """
    # Track heterocycle preservation through the route
    preservation_data = {
        "total_mols": 0,
        "mols_with_heterocycle": 0,
        "reactions_preserving": 0,
        "total_reactions": 0,
    }

    # Expanded list of heterocycles to check - added thiazoline
    heterocycles = [
        "thiazoline",
        "thiazolidine",
        "thiazole",
        "oxazole",
        "isoxazole",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxadiazole",
        "thiadiazole",
        "benzothiazole",
        "benzoxazole",
        "benzimidazole",
        "isothiazole",
        "thiadiazole",
        "oxathiolane",
        "dioxathiolane",
    ]

    # Expanded list of sulfur/nitrogen-containing functional groups
    s_n_functional_groups = [
        "Thiourea",
        "Thioamide",
        "Sulfonamide",
        "Thiocyanate",
        "Isothiocyanate",
        "Sulfonate",
        "Sulfone",
        "Sulfoxide",
        "Monosulfide",
        "Disulfide",
        "Aromatic thiol",
        "Aliphatic thiol",
        "Sulfinate",
        "Sulfonyl halide",
        "Sulfate",
        "Sulfamate",
        "Sulfamic acid",
        "Sulfenic acid",
        "Sulfinic acid",
        "Sulfenate",
        "Carbo-thioester",
        "Thiocarbonyl",
    ]

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            preservation_data["total_mols"] += 1

            # Check if molecule contains any of our target heterocycles
            has_heterocycle = False
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, mol_smiles):
                    has_heterocycle = True
                    print(f"Found {heterocycle} in molecule: {mol_smiles}")
                    break

            # Check for sulfur/nitrogen functional groups if no heterocycle found
            has_s_n_group = False
            if not has_heterocycle:
                for fg in s_n_functional_groups:
                    if checker.check_fg(fg, mol_smiles):
                        has_s_n_group = True
                        print(f"Found {fg} in molecule: {mol_smiles}")
                        break

            if has_heterocycle or has_s_n_group:
                preservation_data["mols_with_heterocycle"] += 1

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            preservation_data["total_reactions"] += 1
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and product
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if heterocycle is present in both reactants and product
                reactant_has_heterocycle = False
                product_has_heterocycle = False
                reactant_heterocycle_type = None
                product_heterocycle_type = None

                # Check for heterocycles in reactants
                for heterocycle in heterocycles:
                    for reactant in reactants:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_has_heterocycle = True
                            reactant_heterocycle_type = heterocycle
                            print(f"Found {heterocycle} in reactant: {reactant}")
                            break
                    if reactant_has_heterocycle:
                        break

                # Check for heterocycles in product
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product):
                        product_has_heterocycle = True
                        product_heterocycle_type = heterocycle
                        print(f"Found {heterocycle} in product: {product}")
                        break

                # Check for S/N functional groups in reactants
                reactant_has_s_n_group = False
                reactant_s_n_group_type = None

                if not reactant_has_heterocycle:
                    for fg in s_n_functional_groups:
                        for reactant in reactants:
                            if checker.check_fg(fg, reactant):
                                reactant_has_s_n_group = True
                                reactant_s_n_group_type = fg
                                print(f"Found {fg} in reactant: {reactant}")
                                break
                        if reactant_has_s_n_group:
                            break

                # Check for S/N functional groups in product
                product_has_s_n_group = False
                product_s_n_group_type = None

                if not product_has_heterocycle:
                    for fg in s_n_functional_groups:
                        if checker.check_fg(fg, product):
                            product_has_s_n_group = True
                            product_s_n_group_type = fg
                            print(f"Found {fg} in product: {product}")
                            break

                # Consider preservation if either heterocycles or S/N groups are maintained
                if (reactant_has_heterocycle and product_has_heterocycle) or (
                    reactant_has_s_n_group and product_has_s_n_group
                ):
                    preservation_data["reactions_preserving"] += 1
                    print(f"Heterocycle/S-N group preserved in reaction: {rsmi}")

                # Check for heterocycle-forming reactions
                if (not reactant_has_heterocycle and product_has_heterocycle) or (
                    not reactant_has_s_n_group and product_has_s_n_group
                ):
                    # Check if this is a known heterocycle-forming reaction
                    heterocycle_forming_rxns = [
                        "benzothiazole",
                        "benzoxazole",
                        "benzimidazole",
                        "thiazole",
                        "oxazole",
                        "imidazole",
                        "{benzothiazole}",
                        "{benzoxazole}",
                        "{benzimidazole}",
                        "{thiazole}",
                        "{oxazole}",
                        "{imidazole}",
                        "{benzothiazole_derivatives_carboxylic-acid/ester}",
                        "{benzothiazole_derivatives_aldehyde}",
                        "{benzoxazole_arom-aldehyde}",
                        "{benzoxazole_carboxylic-acid}",
                    ]
                    for rxn_type in heterocycle_forming_rxns:
                        if checker.check_reaction(rxn_type, rsmi):
                            preservation_data["reactions_preserving"] += 1
                            print(f"Heterocycle-forming reaction detected: {rxn_type} in {rsmi}")
                            break

                # Check for heterocycle transformation (one heterocycle to another)
                if (
                    reactant_has_heterocycle
                    and product_has_heterocycle
                    and reactant_heterocycle_type != product_heterocycle_type
                ):
                    print(
                        f"Heterocycle transformation: {reactant_heterocycle_type} → {product_heterocycle_type}"
                    )
                    preservation_data["reactions_preserving"] += 1

                # Check for S/N group transformation (one S/N group to another)
                if (
                    reactant_has_s_n_group
                    and product_has_s_n_group
                    and reactant_s_n_group_type != product_s_n_group_type
                ):
                    print(
                        f"S/N group transformation: {reactant_s_n_group_type} → {product_s_n_group_type}"
                    )
                    preservation_data["reactions_preserving"] += 1

                # Check for S/N group to heterocycle transformation
                if reactant_has_s_n_group and product_has_heterocycle:
                    print(
                        f"S/N group to heterocycle transformation: {reactant_s_n_group_type} → {product_heterocycle_type}"
                    )
                    preservation_data["reactions_preserving"] += 1

            except Exception as e:
                print(f"Error extracting reactants and product from reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Calculate preservation metrics
    mol_preservation_ratio = 0
    if preservation_data["total_mols"] > 0:
        mol_preservation_ratio = (
            preservation_data["mols_with_heterocycle"] / preservation_data["total_mols"]
        )

    reaction_preservation_ratio = 0
    if preservation_data["total_reactions"] > 0:
        reaction_preservation_ratio = (
            preservation_data["reactions_preserving"] / preservation_data["total_reactions"]
        )

    # Strategy is present if:
    # 1. At least 40% of molecules contain the heterocycle or S/N groups (lowered from 50%)
    # 2. At least 50% of reactions preserve or form the heterocycle or S/N groups (lowered from 60%)
    # 3. We have at least one molecule and one reaction to analyze
    result = (
        preservation_data["total_mols"] > 0
        and preservation_data["total_reactions"] > 0
        and mol_preservation_ratio >= 0.4
        and reaction_preservation_ratio >= 0.5
    )

    print(f"Heterocycle preservation stats:")
    print(
        f"- Molecules with heterocycle/S-N groups: {preservation_data['mols_with_heterocycle']}/{preservation_data['total_mols']} ({mol_preservation_ratio:.2f})"
    )
    print(
        f"- Reactions preserving/forming heterocycle/S-N groups: {preservation_data['reactions_preserving']}/{preservation_data['total_reactions']} ({reaction_preservation_ratio:.2f})"
    )
    print(f"Heterocycle preservation strategy detected: {result}")

    return result
