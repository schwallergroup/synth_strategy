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
    This function detects if the synthetic route uses a formylation reaction
    to introduce an aldehyde group.
    """
    formylation_found = False

    def dfs_traverse(node):
        nonlocal formylation_found

        if node["type"] == "reaction" and not formylation_found:
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains an aldehyde group
            product_has_aldehyde = checker.check_fg("Aldehyde", product_smiles)

            if product_has_aldehyde:
                # Check if any reactant has an aldehyde (to exclude simple aldehyde transfers)
                reactants_with_aldehyde = sum(
                    1 for r in reactants_smiles if checker.check_fg("Aldehyde", r)
                )

                # If product has a new aldehyde (not just transferred from reactants)
                if reactants_with_aldehyde == 0:
                    print(f"Product has new aldehyde: {product_smiles}")

                    # Direct check for known formylation reactions
                    if checker.check_reaction("Friedel-Crafts acylation", rsmi):
                        for reactant in reactants_smiles:
                            if "C(=O)H" in reactant or "CHO" in reactant or "C=O" in reactant:
                                print(f"Formylation detected via Friedel-Crafts acylation: {rsmi}")
                                formylation_found = True
                                return

                    # Check for Vilsmeier-Haack formylation
                    dmf_found = False
                    pocl3_found = False
                    aromatic_found = False

                    for reactant in reactants_smiles:
                        # Check for DMF or similar formamides
                        if (
                            "CN(C)C=O" in reactant
                            or "CN(C)CHO" in reactant
                            or "N(C)C=O" in reactant
                        ):
                            dmf_found = True
                            print(f"DMF or formamide found: {reactant}")

                        # Check for POCl3 or similar reagents
                        if (
                            "POCl3" in reactant
                            or "P(=O)(Cl)(Cl)Cl" in reactant
                            or "PCl" in reactant
                        ):
                            pocl3_found = True
                            print(f"POCl3 or similar found: {reactant}")

                        # Check for aromatic substrate
                        if any(
                            checker.check_ring(ring, reactant)
                            for ring in [
                                "benzene",
                                "pyridine",
                                "furan",
                                "thiophene",
                                "pyrrole",
                                "indole",
                                "quinoline",
                                "isoquinoline",
                            ]
                        ):
                            aromatic_found = True
                            print(f"Aromatic substrate found: {reactant}")

                    if (
                        dmf_found
                        and (pocl3_found or any("Cl" in r for r in reactants_smiles))
                        and aromatic_found
                    ) or (
                        aromatic_found
                        and any("CHO" in r or "C(=O)H" in r for r in reactants_smiles)
                    ):
                        print(f"Vilsmeier-Haack or similar formylation detected: {rsmi}")
                        formylation_found = True
                        return

                    # Check for formylation with CO/H2 (reductive carbonylation)
                    co_found = False
                    h2_or_hydride_found = False
                    metal_catalyst_found = False

                    for reactant in reactants_smiles:
                        if (
                            reactant == "C=O"
                            or reactant == "[C]=O"
                            or reactant == "CO"
                            or "[C-]#[O+]" in reactant
                        ):
                            co_found = True
                            print(f"CO found: {reactant}")
                        if (
                            reactant == "[H][H]"
                            or reactant == "H2"
                            or reactant == "[H2]"
                            or "AlH" in reactant
                            or "BH" in reactant
                        ):
                            h2_or_hydride_found = True
                            print(f"H2 or hydride found: {reactant}")
                        if (
                            "Pd" in reactant
                            or "Rh" in reactant
                            or "Ru" in reactant
                            or "Fe" in reactant
                            or "Co" in reactant
                        ):
                            metal_catalyst_found = True
                            print(f"Metal catalyst found: {reactant}")

                    if co_found and (h2_or_hydride_found or metal_catalyst_found):
                        print(f"Formylation with CO/H2 detected: {rsmi}")
                        formylation_found = True
                        return

                    # Check for formylation via oxidation of methyl groups
                    methyl_oxidation = False
                    for reactant in reactants_smiles:
                        # Check if reactant has methyl group
                        if "C" in reactant and not checker.check_fg("Aldehyde", reactant):
                            # Check if oxidizing agents are present
                            oxidizing_agents = [
                                "KMnO4",
                                "MnO2",
                                "CrO3",
                                "SeO2",
                                "[O]",
                                "O=O",
                                "NaOCl",
                                "H2O2",
                            ]
                            if any(
                                agent in "".join(reactants_smiles) for agent in oxidizing_agents
                            ):
                                methyl_oxidation = True
                                print(f"Potential methyl oxidation to aldehyde: {rsmi}")

                    if methyl_oxidation:
                        print(f"Formylation via methyl oxidation detected: {rsmi}")
                        formylation_found = True
                        return

                    # Check for formylation via hydrolysis of formyl derivatives
                    for reactant in reactants_smiles:
                        if (
                            "OC(=O)H" in reactant
                            or "SC(=O)H" in reactant
                            or "N(C(=O)H)" in reactant
                        ):
                            print(f"Formylation via hydrolysis of formyl derivative: {rsmi}")
                            formylation_found = True
                            return

                    # Check for Duff formylation (hexamine-based)
                    if any("C1N2CN3CN1CN(C2)C3" in r for r in reactants_smiles) and aromatic_found:
                        print(f"Duff formylation detected: {rsmi}")
                        formylation_found = True
                        return

                    # Check for Bouveault aldehyde synthesis (from ester to aldehyde)
                    ester_found = False
                    hydride_found = False

                    for reactant in reactants_smiles:
                        if checker.check_fg("Ester", reactant):
                            ester_found = True
                        if (
                            "AlH" in reactant
                            or "BH" in reactant
                            or "H-" in reactant
                            or "LiAlH" in reactant
                        ):
                            hydride_found = True

                    if ester_found and hydride_found:
                        print(f"Bouveault aldehyde synthesis detected: {rsmi}")
                        formylation_found = True
                        return

                    # Check for Gattermann-Koch formylation
                    if (
                        co_found
                        and any("HCl" in r or "Cl" in r for r in reactants_smiles)
                        and any("AlCl3" in r for r in reactants_smiles)
                        and aromatic_found
                    ):
                        print(f"Gattermann-Koch formylation detected: {rsmi}")
                        formylation_found = True
                        return

                    # Check for reduction of carboxylic acid derivatives to aldehydes
                    carboxylic_derivative_found = False
                    for reactant in reactants_smiles:
                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Nitrile", reactant)
                        ):
                            carboxylic_derivative_found = True

                    if carboxylic_derivative_found and h2_or_hydride_found:
                        print(f"Formylation via reduction of carboxylic acid derivative: {rsmi}")
                        formylation_found = True
                        return

                    # Check for any other reaction that introduces an aldehyde
                    if (
                        checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
                        or checker.check_reaction("Oxidation of alcohol to aldehyde", rsmi)
                        or checker.check_reaction("Bouveault aldehyde synthesis", rsmi)
                        or checker.check_reaction(
                            "Oxidation of aldehydes to carboxylic acids", rsmi
                        )
                    ):  # Reverse direction
                        print(f"Other formylation reaction detected: {rsmi}")
                        formylation_found = True
                        return

                    # General fallback - if we've detected a new aldehyde but none of the specific patterns matched
                    # This is likely a formylation reaction we haven't specifically coded for
                    print(f"General formylation detected (new aldehyde formed): {rsmi}")
                    formylation_found = True
                    return

        # Traverse children
        for child in node.get("children", []):
            if not formylation_found:  # Only continue if we haven't found formylation yet
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return formylation_found
