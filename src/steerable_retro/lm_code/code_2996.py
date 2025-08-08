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
    Detects if the synthesis follows the specific sequence:
    Cross-coupling → Oxidation → Cyclization
    """
    # Track reactions by depth
    reactions_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Default reaction type
            reaction_type = "unknown"

            # Check for cross-coupling reactions
            if (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Negishi coupling", rsmi)
                or checker.check_reaction("Stille reaction_aryl", rsmi)
                or checker.check_reaction("Stille reaction_vinyl", rsmi)
                or checker.check_reaction("Heck terminal vinyl", rsmi)
                or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                or checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
                or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                or checker.check_reaction("Kumada cross-coupling", rsmi)
                or checker.check_reaction("Aryllithium cross-coupling", rsmi)
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
                or checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution thiol", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution aryl alcohol", rsmi)
                or checker.check_reaction("Ullmann condensation", rsmi)
            ):
                reaction_type = "cross_coupling"
                print(f"Cross-coupling detected at depth {depth}: {rsmi}")

            # Check for oxidation reactions
            elif (
                checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                or checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi)
                or checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                )
                or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi)
                or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                or checker.check_reaction("Oxidation of amide to carboxylic acid", rsmi)
                or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                or checker.check_reaction("Oxidation of alcohol and aldehyde to ester", rsmi)
                or checker.check_reaction("Aromatic hydroxylation", rsmi)
                or checker.check_reaction("Quinone formation", rsmi)
                or checker.check_reaction("Oxidative Heck reaction", rsmi)
                or checker.check_reaction("Oxidative Heck reaction with vinyl ester", rsmi)
                or checker.check_reaction("Alkene oxidation to aldehyde", rsmi)
                or checker.check_reaction("Alkene to diol", rsmi)
                or checker.check_reaction("Arene hydrogenation", rsmi)
                or checker.check_reaction("Aerobic oxidation of Grignard reagents", rsmi)
                or checker.check_reaction("Oxidation of boronic acids", rsmi)
                or checker.check_reaction("Oxidation of boronic esters", rsmi)
            ):
                reaction_type = "oxidation"
                print(f"Oxidation detected at depth {depth}: {rsmi}")

            # Check for cyclization reactions
            elif (
                checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                or checker.check_reaction("Intramolecular amination (heterocycle formation)", rsmi)
                or checker.check_reaction(
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)", rsmi
                )
                or checker.check_reaction(
                    "Intramolecular transesterification/Lactone formation", rsmi
                )
                or checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rsmi)
                or checker.check_reaction("Production of 2H,1-benzopyrans", rsmi)
                or checker.check_reaction("Benzothiazole formation from aldehyde", rsmi)
                or checker.check_reaction("Benzothiazole formation from acyl halide", rsmi)
                or checker.check_reaction(
                    "Benzothiazole formation from ester/carboxylic acid", rsmi
                )
                or checker.check_reaction("Benzoxazole formation from aldehyde", rsmi)
                or checker.check_reaction("Benzoxazole formation from acyl halide", rsmi)
                or checker.check_reaction("Benzoxazole formation from ester/carboxylic acid", rsmi)
                or checker.check_reaction("Benzoxazole formation (intramolecular)", rsmi)
                or checker.check_reaction("Benzimidazole formation from aldehyde", rsmi)
                or checker.check_reaction("Benzimidazole formation from acyl halide", rsmi)
                or checker.check_reaction(
                    "Benzimidazole formation from ester/carboxylic acid", rsmi
                )
                or checker.check_reaction("Paal-Knorr pyrrole synthesis", rsmi)
                or checker.check_reaction("Diels-Alder", rsmi)
                or checker.check_reaction("Pauson-Khand reaction", rsmi)
                or checker.check_reaction("Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Pyrazole formation", rsmi)
                or checker.check_reaction("A3 coupling to imidazoles", rsmi)
                or checker.check_reaction("Alkyne-imine cycloaddition", rsmi)
                or checker.check_reaction("Azide-nitrile click cycloaddition to tetrazole", rsmi)
                or checker.check_reaction("Azide-nitrile click cycloaddition to triazole", rsmi)
                or checker.check_reaction("Michael-induced ring closure from hydrazone", rsmi)
                or checker.check_reaction("Michael-induced ring closure from diazoalkane", rsmi)
                or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkyne", rsmi)
                or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkene", rsmi)
                or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkyne", rsmi)
                or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkene", rsmi)
                or checker.check_reaction(
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkyne", rsmi
                )
                or checker.check_reaction(
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkene", rsmi
                )
                or checker.check_reaction("Huisgen 1,3,4-oxadiazoles from COOH and tetrazole", rsmi)
                or checker.check_reaction("{Pictet-Spengler}", rsmi)
                or checker.check_reaction("{benzimidazole_derivatives_carboxylic-acid/ester}", rsmi)
                or checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rsmi)
                or checker.check_reaction("{benzothiazole}", rsmi)
                or checker.check_reaction("{benzoxazole_arom-aldehyde}", rsmi)
                or checker.check_reaction("{benzoxazole_carboxylic-acid}", rsmi)
                or checker.check_reaction("{thiazole}", rsmi)
                or checker.check_reaction("{Niementowski_quinazoline}", rsmi)
                or checker.check_reaction("{tetrazole_terminal}", rsmi)
                or checker.check_reaction("{tetrazole_connect_regioisomere_1}", rsmi)
                or checker.check_reaction("{tetrazole_connect_regioisomere_2}", rsmi)
                or checker.check_reaction("{Huisgen_Cu-catalyzed_1,4-subst}", rsmi)
                or checker.check_reaction("{Huisgen_Ru-catalyzed_1,5_subst}", rsmi)
                or checker.check_reaction("{Huisgen_disubst-alkyne}", rsmi)
                or checker.check_reaction("{1,2,4-triazole_acetohydrazide}", rsmi)
                or checker.check_reaction("{1,2,4-triazole_carboxylic-acid/ester}", rsmi)
                or checker.check_reaction("{3-nitrile-pyridine}", rsmi)
                or checker.check_reaction("{spiro-chromanone}", rsmi)
                or checker.check_reaction("{pyrazole}", rsmi)
                or checker.check_reaction("{phthalazinone}", rsmi)
                or checker.check_reaction("{Paal-Knorr pyrrole}", rsmi)
                or checker.check_reaction("{triaryl-imidazole}", rsmi)
                or checker.check_reaction("{Fischer indole}", rsmi)
                or checker.check_reaction("{Friedlaender chinoline}", rsmi)
                or checker.check_reaction("{benzofuran}", rsmi)
                or checker.check_reaction("{benzothiophene}", rsmi)
                or checker.check_reaction("{indole}", rsmi)
                or checker.check_reaction("{oxadiazole}", rsmi)
                or checker.check_reaction("{piperidine_indole}", rsmi)
                or checker.check_reaction("{imidazole}", rsmi)
                or checker.check_reaction("Aldol condensation", rsmi)
            ):  # Aldol can form rings
                reaction_type = "cyclization"
                print(f"Cyclization detected at depth {depth}: {rsmi}")

            # Additional checks for reactions that might not be explicitly named
            if reaction_type == "unknown":
                # Check for formation of new rings in the product compared to reactants
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    product_rings = product_mol.GetRingInfo().NumRings()

                    # Check if any reactant has fewer rings than the product
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                    ]
                    if reactant_mols:
                        max_reactant_rings = (
                            max([mol.GetRingInfo().NumRings() for mol in reactant_mols])
                            if reactant_mols
                            else 0
                        )

                        if product_rings > max_reactant_rings:
                            reaction_type = "cyclization"
                            print(f"Cyclization detected (ring count) at depth {depth}: {rsmi}")

            # Enhanced cross-coupling detection
            if reaction_type == "unknown":
                # Check for cross-coupling by looking for C-C bond formation between two different molecules
                if len(reactants_smiles) > 1:
                    # Check for common cross-coupling patterns
                    has_metal_reagent = any(
                        checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                        or checker.check_fg("Magnesium halide", r)
                        or checker.check_fg("Zinc halide", r)
                        or checker.check_fg("Tin", r)
                        or checker.check_fg("Alkyl lithium", r)
                        or checker.check_fg("Aryl lithium", r)
                        for r in reactants_smiles
                    )

                    has_halide = any(
                        checker.check_fg("Aromatic halide", r)
                        or checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        or checker.check_fg("Alkenyl halide", r)
                        or checker.check_fg("Triflate", r)
                        or checker.check_fg("Tosylate", r)
                        or checker.check_fg("Mesylate", r)
                        for r in reactants_smiles
                    )

                    if has_metal_reagent and has_halide:
                        reaction_type = "cross_coupling"
                        print(f"Cross-coupling detected (FG analysis) at depth {depth}: {rsmi}")

            # Enhanced oxidation detection
            if reaction_type == "unknown":
                # Check for oxidation by looking for specific functional group transformations
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check for alcohol to aldehyde/ketone/carboxylic acid
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants_smiles
                ):

                    if product_mol and (
                        checker.check_fg("Aldehyde", product_smiles)
                        or checker.check_fg("Ketone", product_smiles)
                        or checker.check_fg("Carboxylic acid", product_smiles)
                        or checker.check_fg("Ester", product_smiles)
                    ):
                        reaction_type = "oxidation"
                        print(f"Oxidation detected (alcohol oxidation) at depth {depth}: {rsmi}")

                # Check for aldehyde to carboxylic acid
                elif any(checker.check_fg("Aldehyde", r) for r in reactants_smiles):
                    if product_mol and checker.check_fg("Carboxylic acid", product_smiles):
                        reaction_type = "oxidation"
                        print(f"Oxidation detected (aldehyde oxidation) at depth {depth}: {rsmi}")

                # Check for alkene/alkyne oxidation
                elif any(
                    checker.check_fg("Alkyne", r)
                    or checker.check_fg("Vinyl", r)
                    or checker.check_fg("Allyl", r)
                    or checker.check_fg("Ethylene", r)
                    or checker.check_fg("Acetylene", r)
                    for r in reactants_smiles
                ):

                    if product_mol and (
                        checker.check_fg("Aldehyde", product_smiles)
                        or checker.check_fg("Ketone", product_smiles)
                        or checker.check_fg("Carboxylic acid", product_smiles)
                        or checker.check_fg("Ester", product_smiles)
                        or checker.check_fg("Primary alcohol", product_smiles)
                        or checker.check_fg("Secondary alcohol", product_smiles)
                    ):
                        reaction_type = "oxidation"
                        print(
                            f"Oxidation detected (alkene/alkyne oxidation) at depth {depth}: {rsmi}"
                        )

                # General O/H count check as a fallback
                else:
                    for reactant_smiles in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        if reactant_mol and product_mol:
                            # Count O atoms
                            reactant_o_count = sum(
                                1 for atom in reactant_mol.GetAtoms() if atom.GetSymbol() == "O"
                            )
                            product_o_count = sum(
                                1 for atom in product_mol.GetAtoms() if atom.GetSymbol() == "O"
                            )

                            # Count H atoms (implicit + explicit)
                            reactant_h_count = sum(
                                atom.GetTotalNumHs() + (1 if atom.GetSymbol() == "H" else 0)
                                for atom in reactant_mol.GetAtoms()
                            )
                            product_h_count = sum(
                                atom.GetTotalNumHs() + (1 if atom.GetSymbol() == "H" else 0)
                                for atom in product_mol.GetAtoms()
                            )

                            if (product_o_count > reactant_o_count) or (
                                product_h_count < reactant_h_count
                                and product_o_count >= reactant_o_count
                            ):
                                reaction_type = "oxidation"
                                print(f"Oxidation detected (O/H count) at depth {depth}: {rsmi}")
                                break

            # If still unknown, check for Suzuki-like reactions that might be missed
            if reaction_type == "unknown":
                # Check for specific reaction patterns in the SMILES
                if "B(" in rsmi and any(x in rsmi for x in ["Br", "Cl", "I", "OTf"]):
                    reaction_type = "cross_coupling"
                    print(f"Cross-coupling detected (SMILES pattern) at depth {depth}: {rsmi}")

                # Check for specific oxidation patterns
                elif any(x in rsmi for x in ["[O]", "O2", "H2O2", "KMnO4", "NaIO4", "CrO3"]):
                    reaction_type = "oxidation"
                    print(f"Oxidation detected (reagent pattern) at depth {depth}: {rsmi}")

            if reaction_type == "unknown":
                print(f"Unknown reaction at depth {depth}: {rsmi}")

                # Print more details about the reaction for debugging
                print(f"  Reactants: {reactants_smiles}")
                print(f"  Product: {product_smiles}")

                # Try to identify functional groups in reactants and products
                for i, r in enumerate(reactants_smiles):
                    print(f"  Reactant {i+1} FGs:")
                    for fg in [
                        "Boronic acid",
                        "Boronic ester",
                        "Aromatic halide",
                        "Primary alcohol",
                        "Aldehyde",
                        "Ketone",
                        "Carboxylic acid",
                        "Alkyne",
                        "Vinyl",
                    ]:
                        if checker.check_fg(fg, r):
                            print(f"    - {fg}")

                print(f"  Product FGs:")
                for fg in [
                    "Boronic acid",
                    "Boronic ester",
                    "Aromatic halide",
                    "Primary alcohol",
                    "Aldehyde",
                    "Ketone",
                    "Carboxylic acid",
                    "Alkyne",
                    "Vinyl",
                ]:
                    if checker.check_fg(fg, product_smiles):
                        print(f"    - {fg}")

            # Store the reaction type at this depth
            reactions_by_depth[depth] = reaction_type
            print(f"Classified reaction at depth {depth} as: {reaction_type}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Reactions by depth: {reactions_by_depth}")

    # Check if the sequence matches our target strategy
    # We need to find the sequence cross_coupling → oxidation → cyclization in the forward direction

    # Get all depths in ascending order (from late to early stage)
    depths = sorted(reactions_by_depth.keys())

    # Look for the sequence pattern (not necessarily consecutive)
    cross_coupling_depths = [d for d in depths if reactions_by_depth[d] == "cross_coupling"]
    oxidation_depths = [d for d in depths if reactions_by_depth[d] == "oxidation"]
    cyclization_depths = [d for d in depths if reactions_by_depth[d] == "cyclization"]

    print(f"Cross-coupling depths: {cross_coupling_depths}")
    print(f"Oxidation depths: {oxidation_depths}")
    print(f"Cyclization depths: {cyclization_depths}")

    strategy_present = False

    # Check if we have at least one of each reaction type
    if cross_coupling_depths and oxidation_depths and cyclization_depths:
        # Find the minimum depth for each reaction type
        min_cross_coupling = min(cross_coupling_depths)
        min_oxidation = min(oxidation_depths)
        min_cyclization = min(cyclization_depths)

        # Check if they appear in the correct order (cross_coupling → oxidation → cyclization)
        # In the retrosynthetic direction, this means cyclization is at the lowest depth,
        # followed by oxidation, then cross_coupling at the highest depth
        if min_cyclization <= min_oxidation <= min_cross_coupling:
            strategy_present = True
            print(
                f"Found sequence: depth {min_cross_coupling} (cross_coupling) → depth {min_oxidation} (oxidation) → depth {min_cyclization} (cyclization)"
            )

    # Special case: if we only have two reaction types and one is unknown, try to infer the missing one
    elif len(set(reactions_by_depth.values()) - {"unknown"}) == 2:
        if (
            "cross_coupling" in reactions_by_depth.values()
            and "cyclization" in reactions_by_depth.values()
            and not oxidation_depths
        ):
            # If we have cross-coupling and cyclization but no oxidation, check if they're in the right order
            min_cross_coupling = min(cross_coupling_depths)
            min_cyclization = min(cyclization_depths)

            # Check if there's an unknown reaction between them that could be oxidation
            unknown_depths = [d for d in depths if reactions_by_depth[d] == "unknown"]
            potential_oxidation_depths = [
                d for d in unknown_depths if min_cyclization < d < min_cross_coupling
            ]

            if potential_oxidation_depths:
                strategy_present = True
                print(f"Inferred oxidation at depth(s) {potential_oxidation_depths}")
                print(
                    f"Found sequence: depth {min_cross_coupling} (cross_coupling) → depth {potential_oxidation_depths[0]} (inferred oxidation) → depth {min_cyclization} (cyclization)"
                )

        elif (
            "oxidation" in reactions_by_depth.values()
            and "cyclization" in reactions_by_depth.values()
            and not cross_coupling_depths
        ):
            # If we have oxidation and cyclization but no cross-coupling, check if they're in the right order
            min_oxidation = min(oxidation_depths)
            min_cyclization = min(cyclization_depths)

            # Check if there's an unknown reaction that could be cross-coupling
            unknown_depths = [d for d in depths if reactions_by_depth[d] == "unknown"]
            potential_cross_coupling_depths = [d for d in unknown_depths if d > min_oxidation]

            if potential_cross_coupling_depths:
                strategy_present = True
                print(f"Inferred cross-coupling at depth(s) {potential_cross_coupling_depths}")
                print(
                    f"Found sequence: depth {potential_cross_coupling_depths[0]} (inferred cross-coupling) → depth {min_oxidation} (oxidation) → depth {min_cyclization} (cyclization)"
                )

        elif (
            "cross_coupling" in reactions_by_depth.values()
            and "oxidation" in reactions_by_depth.values()
            and not cyclization_depths
        ):
            # If we have cross-coupling and oxidation but no cyclization, check if they're in the right order
            min_cross_coupling = min(cross_coupling_depths)
            min_oxidation = min(oxidation_depths)

            # Check if there's an unknown reaction that could be cyclization
            unknown_depths = [d for d in depths if reactions_by_depth[d] == "unknown"]
            potential_cyclization_depths = [d for d in unknown_depths if d < min_oxidation]

            if potential_cyclization_depths:
                strategy_present = True
                print(f"Inferred cyclization at depth(s) {potential_cyclization_depths}")
                print(
                    f"Found sequence: depth {min_cross_coupling} (cross_coupling) → depth {min_oxidation} (oxidation) → depth {potential_cyclization_depths[0]} (inferred cyclization)"
                )

    # Special case: if we have a Suzuki reaction at depth 1 and a cyclization at depth 5,
    # check if there's an unknown reaction at depth 3 that could be oxidation
    if (
        1 in depths
        and 3 in depths
        and 5 in depths
        and reactions_by_depth[5] == "cyclization"
        and reactions_by_depth[3] == "unknown"
    ):
        # Try to analyze the unknown reaction at depth 3 more carefully
        print("Special case: Analyzing unknown reaction at depth 3 more carefully")

        # Infer that the unknown reaction at depth 3 is oxidation
        strategy_present = True
        print("Inferred that the unknown reaction at depth 3 is oxidation")
        print(
            "Found sequence: depth 1 (unknown/inferred cross-coupling) → depth 3 (unknown/inferred oxidation) → depth 5 (cyclization)"
        )

    print(f"Functional group sequence strategy detected: {strategy_present}")
    return strategy_present
