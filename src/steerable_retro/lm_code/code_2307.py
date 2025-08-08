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
    This function detects a synthetic strategy where a benzimidazole ring is formed
    in the final step of the synthesis, preceded by nitro reduction and early biaryl formation.
    """
    # Initialize flags for key features
    benzimidazole_formation_at_depth_1 = False
    nitro_reduction = False
    suzuki_coupling = False
    snar_reaction = False

    # Track molecules through the synthesis to ensure continuity
    molecule_track = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal benzimidazole_formation_at_depth_1, nitro_reduction, suzuki_coupling, snar_reaction, max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Store the relationship between product and reactants
            if product_smiles not in molecule_track:
                molecule_track[product_smiles] = {"depth": depth, "reactants": reactants_smiles}

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            if depth == 1:
                # Check for benzimidazole formation in the final step
                if checker.check_ring("benzimidazole", product_smiles):
                    # Verify it's a formation by checking reactants don't have benzimidazole
                    reactants_have_benzimidazole = any(
                        checker.check_ring("benzimidazole", r) for r in reactants_smiles
                    )

                    # Check for specific benzimidazole formation reactions
                    benzimidazole_reactions = [
                        "benzimidazole_derivatives_aldehyde",
                        "benzimidazole_derivatives_carboxylic-acid/ester",
                        "benzimidazole_formation from aldehyde",
                        "benzimidazole_formation from acyl halide",
                        "benzimidazole_formation from ester/carboxylic acid",
                        "benzimidazole_aldehyde",
                    ]

                    is_benzimidazole_formation = any(
                        checker.check_reaction(rxn, rsmi) for rxn in benzimidazole_reactions
                    )

                    # If no specific reaction is detected, check for general pattern
                    if not is_benzimidazole_formation:
                        # Check if reactants have o-phenylenediamine and aldehyde/carboxylic acid
                        has_diamine = any(
                            checker.check_fg("Primary amine", r) or checker.check_fg("Aniline", r)
                            for r in reactants_smiles
                        )

                        has_carbonyl = any(
                            checker.check_fg("Aldehyde", r)
                            or checker.check_fg("Carboxylic acid", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Acyl halide", r)
                            for r in reactants_smiles
                        )

                        if has_diamine and has_carbonyl:
                            is_benzimidazole_formation = True

                    if not reactants_have_benzimidazole and is_benzimidazole_formation:
                        print("Detected benzimidazole formation at depth 1")
                        benzimidazole_formation_at_depth_1 = True

            else:
                # For reactions at other depths, check if the product is a reactant in a reaction at depth-1
                # This ensures continuity in the synthetic pathway
                is_part_of_pathway = False
                for mol, info in molecule_track.items():
                    if info["depth"] == depth - 1 and product_smiles in info["reactants"]:
                        is_part_of_pathway = True
                        break

                if is_part_of_pathway or depth == max_depth:
                    if depth == 3:
                        # Check for nitro reduction
                        if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                            print("Detected nitro reduction at depth 3")
                            nitro_reduction = True

                        # Alternative check using functional groups if reaction check fails
                        if not nitro_reduction:
                            nitro_in_reactants = any(
                                checker.check_fg("Nitro group", r) for r in reactants_smiles
                            )
                            amine_in_product = checker.check_fg(
                                "Primary amine", product_smiles
                            ) or checker.check_fg("Aniline", product_smiles)

                            if nitro_in_reactants and amine_in_product:
                                print(
                                    "Detected nitro reduction at depth 3 (using functional groups)"
                                )
                                nitro_reduction = True

                    elif depth == 5:
                        # Check for SNAr reaction
                        snar_reactions = [
                            "nucl_sub_aromatic_ortho_nitro",
                            "nucl_sub_aromatic_para_nitro",
                            "heteroaromatic_nuc_sub",
                            "N-arylation",
                            "Ullmann-Goldberg Substitution amine",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                            "Buchwald-Hartwig",
                        ]

                        for rxn_type in snar_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Detected SNAr reaction ({rxn_type}) at depth 5")
                                snar_reaction = True
                                break

                        # Alternative check using functional groups if reaction check fails
                        if not snar_reaction:
                            halide_reactant = any(
                                checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                            )
                            fluoro_reactant = any("F" in r for r in reactants_smiles)
                            amine_reactant = any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Aniline", r)
                                for r in reactants_smiles
                            )

                            # Check if product has a new C-N bond where a C-F or C-halide was
                            if (halide_reactant or fluoro_reactant) and amine_reactant:
                                print("Detected SNAr reaction at depth 5 (using functional groups)")
                                snar_reaction = True

                    elif depth >= 3:
                        # Check for Suzuki coupling
                        suzuki_reactions = [
                            "Suzuki",
                            "Suzuki coupling with boronic acids",
                            "Suzuki coupling with boronic esters",
                            "Suzuki coupling with boronic acids OTf",
                            "Suzuki coupling with boronic esters OTf",
                            "Suzuki coupling with sulfonic esters",
                        ]

                        for rxn_type in suzuki_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Detected Suzuki coupling ({rxn_type}) at depth {depth}")
                                suzuki_coupling = True
                                break

                        # Alternative check using functional groups if reaction check fails
                        if not suzuki_coupling:
                            boronic_acid = any(
                                checker.check_fg("Boronic acid", r) for r in reactants_smiles
                            )
                            boronic_ester = any(
                                checker.check_fg("Boronic ester", r) for r in reactants_smiles
                            )
                            aryl_halide = any(
                                checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                            )

                            if (boronic_acid or boronic_ester) and aryl_halide:
                                print(
                                    f"Detected Suzuki coupling at depth {depth} (using functional groups)"
                                )
                                suzuki_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        benzimidazole_formation_at_depth_1
        and nitro_reduction
        and (suzuki_coupling or snar_reaction)  # Either Suzuki or SNAr should be present
    )

    print(f"Late-stage benzimidazole formation strategy detected: {strategy_present}")
    print(f"Benzimidazole formation: {benzimidazole_formation_at_depth_1}")
    print(f"Nitro reduction: {nitro_reduction}")
    print(f"SNAr reaction: {snar_reaction}")
    print(f"Suzuki coupling: {suzuki_coupling}")

    return strategy_present
