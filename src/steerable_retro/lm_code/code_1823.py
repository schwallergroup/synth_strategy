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
    This function detects if the synthesis involves multiple C-N bond formations.
    Returns True if 3 or more C-N bond formations are detected.
    """
    c_n_bond_formations = 0

    # List of reaction types that typically form C-N bonds
    c_n_bond_forming_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Acylation of secondary amines with anhydrides",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with primary amine to imide",
        "Acyl chloride with secondary amine to amide",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Urea synthesis via isocyanate and diazo",
        "Urea synthesis via isocyanate and sulfonamide",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Ring opening of epoxide with amine",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Aminolysis of esters",
        "Chan-Lam amine",
        "Displacement of ethoxy group by primary amine",
        "Displacement of ethoxy group by secondary amine",
        "aza-Michael addition aromatic",
        "aza-Michael addition secondary",
        "aza-Michael addition primary",
        "Paal-Knorr pyrrole synthesis",
        "Pictet-Spengler",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "pyrazole",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "indole",
        "oxadiazole",
        "imidazole",
        "urea",
        "thiourea",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Mitsunobu_imide",
        "Mitsunobu_sulfonamide",
        "Mitsunobu_tetrazole_1",
        "Mitsunobu_tetrazole_2",
        "Mitsunobu_tetrazole_3",
        "Mitsunobu_tetrazole_4",
        "Eschweiler-Clarke Primary Amine Methylation",
        "Eschweiler-Clarke Secondary Amine Methylation",
        "Reductive methylation of primary amine with formaldehyde",
        "N-methylation",
        "Ynamide synthesis",
        "Phthalimide deprotection",
        "Carboxylic acid to amide conversion",
        "Amine and thiophosgene to isothiocyanate",
        "Phthalic anhydride to phthalimide",
        "Hydrazine synthesis from amine",
    ]

    # Additional heterocyclic rings that involve C-N bond formation
    c_n_heterocycles = [
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
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "acridine",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "pteridin",
        "indazole",
        "benzotriazole",
        "pyrroline",
        "pyrrolidone",
        "imidazolidine",
        "thiazolidine",
        "oxazolidine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # Functional groups that can participate in C-N bond formation
    c_n_reactant_fgs = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Aniline",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Azide",
        "Nitrile",
        "Nitro group",
        "Isocyanate",
        "Isothiocyanate",
        "Hydrazine",
        "Hydrazone",
        "Acylhydrazine",
        "Hydrazone amide",
    ]

    # Functional groups that can be electrophiles in C-N bond formation
    electrophile_fgs = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Acyl halide",
        "Sulfonyl halide",
        "Ester",
        "Carboxylic acid",
        "Anhydride",
        "Aldehyde",
        "Ketone",
        "Epoxide",
        "Aziridine",
        "Oxirane",
        "Triflate",
        "Mesylate",
        "Tosylate",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal c_n_bond_formations

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if the reaction is a known C-N bond forming reaction
                for reaction_type in c_n_bond_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        c_n_bond_formations += 1
                        print(f"C-N bond formation detected: {reaction_type}")
                        break
                else:  # Only execute if no break occurred in the for loop
                    # If no specific reaction type was found, check for general patterns
                    # that might indicate C-N bond formation
                    if ">" in rsmi:
                        reactants = rsmi.split(">")[0]
                        product = rsmi.split(">")[-1]

                        # Check for heterocycle formation (which often involves C-N bond formation)
                        product_has_new_heterocycle = False
                        for ring in c_n_heterocycles:
                            if not checker.check_ring(ring, reactants) and checker.check_ring(
                                ring, product
                            ):
                                c_n_bond_formations += 1
                                print(f"C-N bond formation detected: Formation of {ring} ring")
                                product_has_new_heterocycle = True
                                break

                        if not product_has_new_heterocycle:
                            # Check for Buchwald-Hartwig coupling (general case)
                            if (
                                checker.check_fg("Aromatic halide", reactants)
                                and (
                                    checker.check_fg("Primary amine", reactants)
                                    or checker.check_fg("Secondary amine", reactants)
                                    or checker.check_fg("Aniline", reactants)
                                )
                                and not checker.check_fg("Aromatic halide", product)
                            ):
                                c_n_bond_formations += 1
                                print(
                                    f"C-N bond formation detected: Buchwald-Hartwig coupling (general)"
                                )

                            # Check for reductive amination (general case)
                            elif (
                                (
                                    checker.check_fg("Aldehyde", reactants)
                                    or checker.check_fg("Ketone", reactants)
                                )
                                and (
                                    checker.check_fg("Primary amine", reactants)
                                    or checker.check_fg("Secondary amine", reactants)
                                    or checker.check_fg("Aniline", reactants)
                                )
                                and (
                                    (
                                        not checker.check_fg("Aldehyde", product)
                                        and checker.check_fg("Aldehyde", reactants)
                                    )
                                    or (
                                        not checker.check_fg("Ketone", product)
                                        and checker.check_fg("Ketone", reactants)
                                    )
                                )
                                and (
                                    checker.check_fg("Secondary amine", product)
                                    or checker.check_fg("Tertiary amine", product)
                                )
                            ):
                                c_n_bond_formations += 1
                                print(f"C-N bond formation detected: Reductive amination (general)")

                            # Check for amide formation (general case)
                            elif (
                                (
                                    checker.check_fg("Carboxylic acid", reactants)
                                    or checker.check_fg("Acyl halide", reactants)
                                    or checker.check_fg("Ester", reactants)
                                    or checker.check_fg("Anhydride", reactants)
                                )
                                and (
                                    checker.check_fg("Primary amine", reactants)
                                    or checker.check_fg("Secondary amine", reactants)
                                    or checker.check_fg("Aniline", reactants)
                                )
                                and (
                                    checker.check_fg("Primary amide", product)
                                    or checker.check_fg("Secondary amide", product)
                                    or checker.check_fg("Tertiary amide", product)
                                )
                            ):
                                c_n_bond_formations += 1
                                print(f"C-N bond formation detected: Amide formation (general)")

                            # Check for urea formation (general case)
                            elif (
                                checker.check_fg("Isocyanate", reactants)
                                and (
                                    checker.check_fg("Primary amine", reactants)
                                    or checker.check_fg("Secondary amine", reactants)
                                    or checker.check_fg("Aniline", reactants)
                                )
                                and checker.check_fg("Urea", product)
                            ):
                                c_n_bond_formations += 1
                                print(f"C-N bond formation detected: Urea formation (general)")

                            # Check for sulfonamide formation (general case)
                            elif (
                                checker.check_fg("Sulfonyl halide", reactants)
                                and (
                                    checker.check_fg("Primary amine", reactants)
                                    or checker.check_fg("Secondary amine", reactants)
                                    or checker.check_fg("Aniline", reactants)
                                )
                                and checker.check_fg("Sulfonamide", product)
                            ):
                                c_n_bond_formations += 1
                                print(
                                    f"C-N bond formation detected: Sulfonamide formation (general)"
                                )

                            # Check for imine formation (general case)
                            elif (
                                (
                                    checker.check_fg("Aldehyde", reactants)
                                    or checker.check_fg("Ketone", reactants)
                                )
                                and checker.check_fg("Primary amine", reactants)
                                and (
                                    checker.check_fg("Substituted imine", product)
                                    or checker.check_fg("Unsubstituted imine", product)
                                )
                            ):
                                c_n_bond_formations += 1
                                print(f"C-N bond formation detected: Imine formation (general)")

                            # Check for aza-Michael addition (general case)
                            elif (
                                checker.check_fg("Primary amine", reactants)
                                or checker.check_fg("Secondary amine", reactants)
                                or checker.check_fg("Aniline", reactants)
                            ) and (
                                checker.check_fg("Secondary amine", product)
                                or checker.check_fg("Tertiary amine", product)
                            ):
                                # Check if there's a conjugated system in reactants
                                for fg in ["Vinyl", "Allyl", "Conjugated diene"]:
                                    if checker.check_fg(fg, reactants) and not checker.check_fg(
                                        fg, product
                                    ):
                                        c_n_bond_formations += 1
                                        print(
                                            f"C-N bond formation detected: aza-Michael addition (general)"
                                        )
                                        break

                            # Check for alkylation of amines (general case)
                            elif (
                                (
                                    checker.check_fg("Primary halide", reactants)
                                    or checker.check_fg("Secondary halide", reactants)
                                    or checker.check_fg("Tertiary halide", reactants)
                                    or checker.check_fg("Triflate", reactants)
                                    or checker.check_fg("Mesylate", reactants)
                                    or checker.check_fg("Tosylate", reactants)
                                )
                                and (
                                    checker.check_fg("Primary amine", reactants)
                                    or checker.check_fg("Secondary amine", reactants)
                                    or checker.check_fg("Aniline", reactants)
                                )
                                and (
                                    checker.check_fg("Secondary amine", product)
                                    or checker.check_fg("Tertiary amine", product)
                                )
                            ):
                                c_n_bond_formations += 1
                                print(
                                    f"C-N bond formation detected: Alkylation of amines (general)"
                                )

                            # Check for azide to amine conversion (often involves C-N bond formation)
                            elif (
                                checker.check_fg("Azide", reactants)
                                and (
                                    checker.check_fg("Primary amine", product)
                                    or checker.check_fg("Secondary amine", product)
                                    or checker.check_fg("Tertiary amine", product)
                                )
                                and not checker.check_fg("Azide", product)
                            ):
                                c_n_bond_formations += 1
                                print(f"C-N bond formation detected: Azide to amine conversion")

                            # Check for nitrile to amide/amine conversion
                            elif (
                                checker.check_fg("Nitrile", reactants)
                                and (
                                    checker.check_fg("Primary amide", product)
                                    or checker.check_fg("Secondary amide", product)
                                    or checker.check_fg("Primary amine", product)
                                    or checker.check_fg("Secondary amine", product)
                                )
                                and not checker.check_fg("Nitrile", product)
                            ):
                                c_n_bond_formations += 1
                                print(
                                    f"C-N bond formation detected: Nitrile to amide/amine conversion"
                                )

                            # Check for nitro to amine reduction (can lead to C-N bond formation)
                            elif (
                                checker.check_fg("Nitro group", reactants)
                                and (
                                    checker.check_fg("Primary amine", product)
                                    or checker.check_fg("Aniline", product)
                                )
                                and not checker.check_fg("Nitro group", product)
                            ):
                                c_n_bond_formations += 1
                                print(f"C-N bond formation detected: Nitro to amine reduction")

                            # General check for new C-N bonds by comparing nitrogen-containing FGs
                            else:
                                # Count nitrogen-containing FGs in reactants and products
                                reactant_n_fgs = []
                                product_n_fgs = []

                                for fg in c_n_reactant_fgs:
                                    if checker.check_fg(fg, reactants):
                                        reactant_n_fgs.append(fg)
                                    if checker.check_fg(fg, product):
                                        product_n_fgs.append(fg)

                                # If there are more nitrogen FGs in product, likely C-N bond formation
                                if len(product_n_fgs) > len(reactant_n_fgs):
                                    c_n_bond_formations += 1
                                    print(
                                        f"C-N bond formation detected: Increase in nitrogen-containing functional groups"
                                    )
                                    print(f"  Reactant FGs: {reactant_n_fgs}")
                                    print(f"  Product FGs: {product_n_fgs}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Total C-N bond formations detected: {c_n_bond_formations}")
    return c_n_bond_formations >= 3  # Return True if 3 or more C-N bond formations
