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


def main(node, depth, data):
    """
    Traverse the synthesis route using depth-first search and collect data about
    chemical transformations and synthesis structure.
    """
    # Update maximum depth
    data["max_depth"] = max(data["max_depth"], depth)

    # Process reaction nodes
    if node["type"] == "reaction":
        try:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for heterocycle formation
            heterocycles = [
                "furan",
                "pyrrole",
                "thiophene",
                "oxazole",
                "thiazole",
                "pyridine",
                "imidazole",
                "pyrazole",
                "isoxazole",
                "isothiazole",
                "oxadiazole",
                "thiadiazole",
                "triazole",
                "tetrazole",
                "benzoxazole",
                "benzothiazole",
                "benzimidazole",
                "indole",
                "quinoline",
                "isoquinoline",
            ]

            product_has_heterocycle = False
            heterocycle_found = None

            for ring in heterocycles:
                if checker.check_ring(ring, product):
                    product_has_heterocycle = True
                    heterocycle_found = ring
                    break

            reactants_have_heterocycle = any(
                any(checker.check_ring(ring, r) for ring in heterocycles) for r in reactants
            )

            # Check for heterocycle formation reactions
            heterocycle_formation_rxn = any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Formation of NOS Heterocycles",
                    "Paal-Knorr pyrrole synthesis",
                    "{benzimidazole_derivatives_carboxylic-acid/ester}",
                    "{benzimidazole_derivatives_aldehyde}",
                    "{benzothiazole}",
                    "{benzoxazole_arom-aldehyde}",
                    "{benzoxazole_carboxylic-acid}",
                    "{thiazole}",
                    "{tetrazole_terminal}",
                    "{tetrazole_connect_regioisomere_1}",
                    "{tetrazole_connect_regioisomere_2}",
                    "{1,2,4-triazole_acetohydrazide}",
                    "{1,2,4-triazole_carboxylic-acid/ester}",
                    "{pyrazole}",
                    "{Paal-Knorr pyrrole}",
                    "{oxadiazole}",
                    "{indole}",
                    "{Fischer indole}",
                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                    "Huisgen 1,3 dipolar cycloaddition",
                    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                    "Pyrazole formation",
                ]
            )

            if (
                product_has_heterocycle and not reactants_have_heterocycle
            ) or heterocycle_formation_rxn:
                data["heterocycle_formation"] = True
                data["heterocycle_formation_depth"] = depth
                data["heterocycle_type"] = heterocycle_found
                print(f"Detected heterocycle formation ({heterocycle_found}) at depth {depth}")

            # Check for oxime reduction
            reduction_rxns = [
                "Reduction of nitrile to amine",
                "Reduction of aldehydes and ketones to alcohols",
                "Reduction of primary amides to amines",
                "Reduction of secondary amides to amines",
                "Reduction of tertiary amides to amines",
                "Reduction of ester to primary alcohol",
                "Reduction of carboxylic acid to primary alcohol",
                "Azide to amine reduction (Staudinger)",
            ]

            is_reduction = any(checker.check_reaction(rxn, rsmi) for rxn in reduction_rxns)

            if is_reduction:
                product_has_oxime = checker.check_fg("Oxime", product)
                reactants_have_oxime = any(checker.check_fg("Oxime", r) for r in reactants)

                if reactants_have_oxime and not product_has_oxime:
                    data["oxime_reduction"] = True
                    data["oxime_reduction_depth"] = depth
                    print(f"Detected oxime reduction at depth {depth}")

            # Check for halide to alcohol conversion
            halide_alcohol_rxns = [
                "Primary halide to alcohol",
                "Secondary halide to alcohol",
                "Alcohol to chloride_Other",
                "Alcohol to chloride_HCl",
                "Alcohol to chloride_SOCl2",
                "Alcohol to chloride_POCl3",
                "Alcohol to chloride_PCl5_ortho",
                "Alcohol to chloride_POCl3_ortho",
                "Alcohol to chloride_POCl3_para",
            ]

            is_halide_alcohol_rxn = any(
                checker.check_reaction(rxn, rsmi) for rxn in halide_alcohol_rxns
            )

            if is_halide_alcohol_rxn:
                product_has_alcohol = (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                    or checker.check_fg("Aromatic alcohol", product)
                )

                reactants_have_halide = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    or checker.check_fg("Aromatic halide", r)
                    or checker.check_fg("Alkenyl halide", r)
                    for r in reactants
                )

                if product_has_alcohol and reactants_have_halide:
                    data["halide_to_alcohol"] = True
                    data["halide_to_alcohol_depth"] = depth
                    print(f"Detected halide to alcohol conversion at depth {depth}")

            # Check for other important functional group transformations
            important_rxns = [
                "Williamson Ether Synthesis",
                "Esterification of Carboxylic Acids",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
                "Reductive amination with alcohol",
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic esters",
                "Heck terminal vinyl",
                "Oxidative Heck reaction",
                "Heck reaction with vinyl ester and amine",
                "Wittig reaction with triphenylphosphorane",
                "Wittig with Phosphonium",
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Alkylation of amines",
                "Oxidation of alcohol to carboxylic acid",
                "Oxidation of aldehydes to carboxylic acids",
                "Acylation of olefines by aldehydes",
                "Acylation of secondary amines with anhydrides",
                "Schotten-Baumann to ester",
                "{Schotten-Baumann_amide}",
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Aldol condensation",
                "Michael addition",
                "aza-Michael addition primary",
                "aza-Michael addition secondary",
                "Grignard from aldehyde to alcohol",
                "Grignard from ketone to alcohol",
            ]

            if any(checker.check_reaction(rxn, rsmi) for rxn in important_rxns):
                data["other_fg_transformation"] = True
                data["other_fg_transformation_depth"] = depth
                print(f"Detected important functional group transformation at depth {depth}")

            # Check for direct functional group transformations by comparing reactants and products
            if not data["other_fg_transformation"]:
                important_fgs = [
                    "Primary alcohol",
                    "Secondary alcohol",
                    "Tertiary alcohol",
                    "Primary amine",
                    "Secondary amine",
                    "Tertiary amine",
                    "Carboxylic acid",
                    "Ester",
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide",
                    "Nitrile",
                    "Aldehyde",
                    "Ketone",
                    "Ether",
                    "Acyl halide",
                    "Anhydride",
                ]

                for fg in important_fgs:
                    product_has_fg = checker.check_fg(fg, product)
                    reactants_have_fg = any(checker.check_fg(fg, r) for r in reactants)

                    if (product_has_fg and not reactants_have_fg) or (
                        not product_has_fg and reactants_have_fg
                    ):
                        data["other_fg_transformation"] = True
                        data["other_fg_transformation_depth"] = depth
                        print(f"Detected {fg} transformation at depth {depth}")
                        break

        except Exception as e:
            print(f"Error processing reaction node: {e}")

    # Count children for branching factor analysis
    children_count = len(node.get("children", []))
    if children_count > 0:
        data["branching_factor"].append(children_count)

    # Recursively process children
    for child in node.get("children", []):
        dfs_traverse(child, depth + 1, data)
