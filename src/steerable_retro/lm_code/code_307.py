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
    This function detects a synthesis strategy involving a late-stage Suzuki coupling
    as the final step in a linear heterocycle elaboration sequence.
    """
    # Track key reactions and features
    suzuki_reaction_found = False
    suzuki_depth = float("inf")
    alkylation_found = False
    alkylation_depth = float("inf")
    chlorination_found = False
    chlorination_depth = float("inf")
    has_heterocycle_in_final = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_reaction_found, suzuki_depth, alkylation_found, alkylation_depth
        nonlocal chlorination_found, chlorination_depth, has_heterocycle_in_final

        # Check for heterocycle in molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for common heterocycles - expanded list
            heterocycle_rings = [
                "pyridine",
                "pyrimidine",
                "pyrazine",
                "pyridazine",
                "triazole",
                "tetrazole",
                "furan",
                "thiophene",
                "pyrrole",
                "imidazole",
                "oxazole",
                "thiazole",
                "indole",
                "benzimidazole",
                "benzoxazole",
                "benzothiazole",
                "quinoline",
                "isoquinoline",
                "piperidine",
                "piperazine",
                "morpholine",
                "thiomorpholine",
                "pyrrolidine",
                "diazepane",
                "oxadiazole",
                "thiadiazole",
                "purine",
                "carbazole",
                "acridine",
            ]
            has_heterocycle = any(
                checker.check_ring(ring, mol_smiles) for ring in heterocycle_rings
            )

            # If this is the final product (depth 0), record heterocycle presence
            if depth == 0 and has_heterocycle:
                has_heterocycle_in_final = True
                print(f"Found heterocycle in final product: {mol_smiles}")

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for Suzuki coupling - expanded list
                suzuki_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with boronic esters",
                    "{Suzuki}",
                ]

                # Check if this is a Suzuki coupling reaction
                is_suzuki = any(checker.check_reaction(rxn, rsmi) for rxn in suzuki_reactions)

                # Additional check for boronic acid/ester and aryl halide functional groups
                has_boronic = any(
                    checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                    for r in reactants_smiles
                )
                has_aryl_halide = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                )

                if is_suzuki or (has_boronic and has_aryl_halide):
                    print(f"Found Suzuki coupling at depth {depth}: {rsmi}")
                    suzuki_reaction_found = True
                    suzuki_depth = min(suzuki_depth, depth)

                # Check for alkylation reactions - expanded list
                alkylation_reactions = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Friedel-Crafts alkylation",
                    "Friedel-Crafts alkylation with halide",
                    "S-alkylation of thiols",
                    "S-alkylation of thiols (ethyl)",
                    "S-alkylation of thiols with alcohols",
                    "S-alkylation of thiols with alcohols (ethyl)",
                    "Williamson Ether Synthesis",
                    "{Williamson ether}",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "O-alkylation of amides with diazo compounds",
                    "Alkylation of amines",
                    "C-alkylation",
                    "N-methylation",
                    "O-methylation",
                    "S-methylation",
                    "Methylation",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "Methylation with MeI_aryl",
                    "Methylation with MeI_SH",
                    "Methylation with DMS",
                    "Methylation of OH with DMS",
                    "Methylation with DMC",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                ]

                # Additional check for alkylation by looking at functional groups
                has_alkyl_halide = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    for r in reactants_smiles
                )
                has_amine = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    for r in reactants_smiles
                )

                if any(checker.check_reaction(rxn, rsmi) for rxn in alkylation_reactions) or (
                    has_alkyl_halide and has_amine
                ):
                    print(f"Found alkylation reaction at depth {depth}: {rsmi}")
                    alkylation_found = True
                    alkylation_depth = min(alkylation_depth, depth)

                # Check for chlorination reactions - expanded list
                chlorination_reactions = [
                    "Aromatic chlorination",
                    "Chlorination",
                    "Alcohol to chloride_POCl3",
                    "Alcohol to chloride_POCl3_ortho",
                    "Alcohol to chloride_POCl3_para",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_HCl",
                    "Alcohol to chloride_Salt",
                    "Alcohol to chloride_Other",
                    "Alcohol to chloride_sulfonyl chloride",
                    "Alcohol to chloride_CHCl3",
                    "Alcohol to chloride_CH2Cl2",
                    "Alcohol to chloride_PCl5_ortho",
                    "Primary amine to chloride",
                    "Aromatic substitution of bromine by chlorine",
                ]

                # Additional check for chlorination by looking at functional groups
                has_chlorine_in_product = "Cl" in product_smiles
                has_alcohol = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants_smiles
                )

                if any(checker.check_reaction(rxn, rsmi) for rxn in chlorination_reactions) or (
                    has_alcohol and has_chlorine_in_product
                ):
                    print(f"Found chlorination reaction at depth {depth}: {rsmi}")
                    chlorination_found = True
                    chlorination_depth = min(chlorination_depth, depth)

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present with correct sequence
    # Suzuki should be at the latest stage (lowest depth)
    # Either alkylation or chlorination should be present in the sequence
    strategy_present = (
        suzuki_reaction_found
        and has_heterocycle_in_final
        and (alkylation_found or chlorination_found)
        and (suzuki_depth <= 2)  # Suzuki should be in the late stage (depth 0, 1, or 2)
    )

    # Additional check for the sequence - more flexible
    if strategy_present and alkylation_found and chlorination_found:
        # If both alkylation and chlorination are found, at least one should be before Suzuki
        strategy_present = alkylation_depth > suzuki_depth or chlorination_depth > suzuki_depth

    if strategy_present:
        print("Detected late-stage Suzuki coupling strategy with heterocycle elaboration")
    else:
        print(
            f"Strategy detection failed: Suzuki={suzuki_reaction_found} (depth={suzuki_depth}), "
            f"Alkylation={alkylation_found} (depth={alkylation_depth}), "
            f"Chlorination={chlorination_found} (depth={chlorination_depth}), "
            f"Heterocycle={has_heterocycle_in_final}"
        )

    return strategy_present
