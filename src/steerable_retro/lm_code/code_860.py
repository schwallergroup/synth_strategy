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
    This function detects if the synthetic route involves early heterocycle formation
    followed by late-stage amide formation.
    """
    heterocycle_formation_depth = -1
    amide_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, amide_formation_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for heterocycle formation reactions
                heterocycle_reactions = [
                    "benzimidazole_derivatives_carboxylic-acid/ester",
                    "benzimidazole_derivatives_aldehyde",
                    "benzoxazole_arom-aldehyde",
                    "benzoxazole_carboxylic-acid",
                    "benzothiazole",
                    "thiazole",
                    "Paal-Knorr pyrrole",
                    "Fischer indole",
                    "pyrazole",
                    "tetrazole_terminal",
                    "tetrazole_connect_regioisomere_1",
                    "tetrazole_connect_regioisomere_2",
                    "1,2,4-triazole_acetohydrazide",
                    "1,2,4-triazole_carboxylic-acid/ester",
                    "oxadiazole",
                    "imidazole",
                    "Niementowski_quinazoline",
                    "Formation of NOS Heterocycles",
                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                    "Huisgen 1,3 dipolar cycloaddition",
                    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                    "Pyrazole formation",
                    "Huisgen_Cu-catalyzed_1,4-subst",
                    "Huisgen_Ru-catalyzed_1,5_subst",
                    "Huisgen_disubst-alkyne",
                ]

                for reaction_name in heterocycle_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(f"Heterocycle formation detected at depth {depth}: {reaction_name}")
                        if heterocycle_formation_depth == -1 or depth > heterocycle_formation_depth:
                            heterocycle_formation_depth = depth
                        break

                # If no specific heterocycle reaction detected, check product for new heterocycle
                if heterocycle_formation_depth == -1:
                    product = rsmi.split(">")[-1]
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if product contains heterocycle
                    heterocycle_rings = [
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                        "pyrrole",
                        "pyridine",
                        "furan",
                        "thiophene",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrimidine",
                        "triazole",
                        "tetrazole",
                        "indole",
                        "pyrazole",
                        "isoxazole",
                        "isothiazole",
                        "oxadiazole",
                        "thiadiazole",
                        "quinoline",
                        "isoquinoline",
                        "purine",
                        "piperidine",
                        "piperazine",
                        "morpholine",
                        "thiomorpholine",
                        "pyrrolidine",
                        "azetidine",
                        "aziridine",
                        "azepane",
                        "diazepane",
                        "carbazole",
                        "acridine",
                        "benzotriazole",
                        "indazole",
                        "pteridin",
                        "phenothiazine",
                        "phenoxazine",
                    ]

                    for ring_name in heterocycle_rings:
                        if checker.check_ring(ring_name, product):
                            # Check if reactants don't have this heterocycle
                            has_ring_in_reactants = False
                            for reactant in reactants:
                                if checker.check_ring(ring_name, reactant):
                                    has_ring_in_reactants = True
                                    break

                            if not has_ring_in_reactants:
                                print(
                                    f"Heterocycle formation detected at depth {depth}: {ring_name} ring"
                                )
                                if (
                                    heterocycle_formation_depth == -1
                                    or depth > heterocycle_formation_depth
                                ):
                                    heterocycle_formation_depth = depth
                                break

                # Check for amide formation reactions
                amide_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Schotten-Baumann_amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Ester with ammonia to amide",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to imide",
                    "Nitrile to amide",
                    "Hydroxamic Synthesis",
                ]

                for reaction_name in amide_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(f"Amide formation detected at depth {depth}: {reaction_name}")
                        if amide_formation_depth == -1 or depth < amide_formation_depth:
                            amide_formation_depth = depth
                        break

                # If no specific amide reaction detected, check for amide formation by functional group change
                if amide_formation_depth == -1:
                    product = rsmi.split(">")[-1]
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if product has amide but reactants don't
                    if (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        has_amide_in_reactants = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Primary amide", reactant)
                                or checker.check_fg("Secondary amide", reactant)
                                or checker.check_fg("Tertiary amide", reactant)
                            ):
                                has_amide_in_reactants = True
                                break

                        if not has_amide_in_reactants:
                            # Check if reactants have amine and carboxylic acid/derivative
                            has_amine = False
                            has_acid_or_derivative = False

                            for reactant in reactants:
                                if (
                                    checker.check_fg("Primary amine", reactant)
                                    or checker.check_fg("Secondary amine", reactant)
                                    or checker.check_fg("Tertiary amine", reactant)
                                    or checker.check_fg("Aniline", reactant)
                                ):
                                    has_amine = True

                                if (
                                    checker.check_fg("Carboxylic acid", reactant)
                                    or checker.check_fg("Acyl halide", reactant)
                                    or checker.check_fg("Ester", reactant)
                                    or checker.check_fg("Anhydride", reactant)
                                    or checker.check_fg("Nitrile", reactant)
                                ):
                                    has_acid_or_derivative = True

                            if has_amine and has_acid_or_derivative:
                                print(f"Amide formation detected at depth {depth} (by FG analysis)")
                                if amide_formation_depth == -1 or depth < amide_formation_depth:
                                    amide_formation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Heterocycle formation depth: {heterocycle_formation_depth}")
    print(f"Amide formation depth: {amide_formation_depth}")

    # In retrosynthesis, lower depth = later stage, higher depth = earlier stage
    # We want heterocycle formation to be early (higher depth) and amide formation to be late (lower depth)
    if (
        heterocycle_formation_depth > amide_formation_depth
        and heterocycle_formation_depth != -1
        and amide_formation_depth != -1
    ):
        print(
            "Strategy detected: Early heterocycle formation followed by late-stage amide formation"
        )
        return True
    else:
        print("Strategy not detected")
        return False
