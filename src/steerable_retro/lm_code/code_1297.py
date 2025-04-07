#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis involves multiple heteroatom bond formations
    (C-O, C-N, C-S, N-S bonds).
    """
    c_o_bond = False
    c_n_bond = False
    c_s_bond = False
    n_s_bond = False

    # Reactions that form C-O bonds (expanded list)
    c_o_reactions = [
        "Williamson Ether Synthesis",
        "Esterification of Carboxylic Acids",
        "Alcohol protection with silyl ethers",
        "Acetal hydrolysis to diol",
        "Aldehyde or ketone acetalization",
        "Oxidative esterification of primary alcohols",
        "Alcohol to ether",
        "Mitsunobu aryl ether",
        "Mitsunobu esterification",
        "Transesterification",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Schotten-Baumann to ester",
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of ketone to carboxylic acid",
        "Reduction of ester to primary alcohol",
        "Reduction of ketone to secondary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Primary alkyl halide to alcohol",
        "Secondary alkyl halide to alcohol",
        "Markovnikov alkene hydration to alcohol",
        "anti-Markovnikov alkene hydration to alcohol",
        "Alkene to diol",
        "Alkene oxidation to aldehyde",
        "Oxidation of alcohol and aldehyde to ester",
        "Acetic anhydride and alcohol to ester",
        "Chan-Lam etherification",
        "Mitsunobu_phenole",
        "oxa-Michael addition",
    ]

    # Reactions that form C-N bonds (expanded list)
    c_n_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Reduction of nitrile to amine",
        "Reduction of nitro groups to amines",
        "Reduction of primary amides to amines",
        "Reduction of secondary amides to amines",
        "Reduction of tertiary amides to amines",
        "Azide to amine reduction (Staudinger)",
        "Boc amine protection",
        "Boc amine deprotection",
        "Phthalimide deprotection",
        "N-glutarimide deprotection",
        "Eschweiler-Clarke Primary Amine Methylation",
        "Eschweiler-Clarke Secondary Amine Methylation",
        "Reductive methylation of primary amine with formaldehyde",
        "N-methylation",
        "aza-Michael addition aromatic",
        "aza-Michael addition secondary",
        "aza-Michael addition primary",
        "Buchwald-Hartwig",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Ullmann-Goldberg Substitution amine",
        "Carboxylic acid to amide conversion",
    ]

    # Reactions that form C-S bonds
    c_s_reactions = [
        "S-alkylation of thiols",
        "S-alkylation of thiols (ethyl)",
        "S-alkylation of thiols with alcohols",
        "S-alkylation of thiols with alcohols (ethyl)",
        "thia-Michael addition",
        "Ullmann-Goldberg Substitution thiol",
        "S-methylation",
        "Thiol-ene reaction",
    ]

    # Reactions that form N-S bonds
    n_s_reactions = [
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Formation of Sulfonic Esters",
        "Schotten-Baumann to ester",
        "sulfon_amide",
        "Sulfamoylarylamides from carboxylic acids and amines",
    ]

    # C-O functional groups
    c_o_fgs = [
        "Ether",
        "Ester",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Aromatic alcohol",
        "Phenol",
        "Carboxylic acid",
        "Anhydride",
        "Acetal/Ketal",
        "Carbonic Ester",
        "Carbamic ester",
    ]

    # C-N functional groups
    c_n_fgs = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Aniline",
        "Urea",
        "Carboxylic acid",
        "Imine",
        "Nitrile",
        "Azide",
        "Isocyanate",
    ]

    # C-S functional groups
    c_s_fgs = [
        "Aromatic thiol",
        "Aliphatic thiol",
        "Monosulfide",
        "Disulfide",
        "Thioamide",
        "Carbo-thioester",
        "Thiocarbonyl",
    ]

    # N-S functional groups
    n_s_fgs = ["Sulfonamide", "Sulfamate", "Sulfamic acid"]

    def dfs_traverse(node):
        nonlocal c_o_bond, c_n_bond, c_s_bond, n_s_bond

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for C-O bond formation
            if not c_o_bond:
                for reaction in c_o_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        print(f"Detected C-O bond formation via {reaction}")
                        c_o_bond = True
                        break

                # Additional check for C-O functional group formation
                if not c_o_bond:
                    for fg in c_o_fgs:
                        if checker.check_fg(fg, product) and not any(
                            checker.check_fg(fg, r) for r in reactants
                        ):
                            print(f"Detected C-O bond formation via {fg} formation")
                            c_o_bond = True
                            break

            # Check for C-N bond formation
            if not c_n_bond:
                for reaction in c_n_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        print(f"Detected C-N bond formation via {reaction}")
                        c_n_bond = True
                        break

                # Additional check for C-N functional group formation
                if not c_n_bond:
                    for fg in c_n_fgs:
                        if checker.check_fg(fg, product) and not any(
                            checker.check_fg(fg, r) for r in reactants
                        ):
                            print(f"Detected C-N bond formation via {fg} formation")
                            c_n_bond = True
                            break

            # Check for C-S bond formation
            if not c_s_bond:
                for reaction in c_s_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        print(f"Detected C-S bond formation via {reaction}")
                        c_s_bond = True
                        break

                # Additional check for C-S functional group formation
                if not c_s_bond:
                    for fg in c_s_fgs:
                        if checker.check_fg(fg, product) and not any(
                            checker.check_fg(fg, r) for r in reactants
                        ):
                            print(f"Detected C-S bond formation via {fg} formation")
                            c_s_bond = True
                            break

            # Check for N-S bond formation
            if not n_s_bond:
                for reaction in n_s_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        print(f"Detected N-S bond formation via {reaction}")
                        n_s_bond = True
                        break

                # Additional check for N-S functional group formation
                if not n_s_bond:
                    for fg in n_s_fgs:
                        if checker.check_fg(fg, product) and not any(
                            checker.check_fg(fg, r) for r in reactants
                        ):
                            print(f"Detected N-S bond formation via {fg} formation")
                            n_s_bond = True
                            break

            # Check for ring formations that might involve heteroatom bonds
            if not (c_o_bond and c_n_bond and c_s_bond and n_s_bond):
                hetero_rings = [
                    "furan",
                    "pyran",
                    "dioxane",
                    "tetrahydrofuran",
                    "tetrahydropyran",
                    "oxirane",
                    "oxetane",
                    "pyrrole",
                    "pyridine",
                    "pyrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "pyrimidine",
                    "triazole",
                    "tetrazole",
                    "thiophene",
                    "thiopyran",
                    "thiirane",
                    "thietane",
                ]

                for ring in hetero_rings:
                    if checker.check_ring(ring, product) and not any(
                        checker.check_ring(ring, r) for r in reactants
                    ):
                        if (
                            ring
                            in [
                                "furan",
                                "pyran",
                                "dioxane",
                                "tetrahydrofuran",
                                "tetrahydropyran",
                                "oxirane",
                                "oxetane",
                                "oxazole",
                            ]
                            and not c_o_bond
                        ):
                            print(f"Detected C-O bond formation via {ring} ring formation")
                            c_o_bond = True
                        elif (
                            ring
                            in [
                                "pyrrole",
                                "pyridine",
                                "pyrazole",
                                "imidazole",
                                "pyrimidine",
                                "triazole",
                                "tetrazole",
                            ]
                            and not c_n_bond
                        ):
                            print(f"Detected C-N bond formation via {ring} ring formation")
                            c_n_bond = True
                        elif (
                            ring in ["thiophene", "thiopyran", "thiirane", "thietane", "thiazole"]
                            and not c_s_bond
                        ):
                            print(f"Detected C-S bond formation via {ring} ring formation")
                            c_s_bond = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if at least two different heteroatom bond formations were detected
    bond_count = sum([c_o_bond, c_n_bond, c_s_bond, n_s_bond])
    result = bond_count >= 2

    print(f"C-O bond: {c_o_bond}, C-N bond: {c_n_bond}, C-S bond: {c_s_bond}, N-S bond: {n_s_bond}")
    print(f"Multiple heteroatom functionalization: {result}")

    return result
