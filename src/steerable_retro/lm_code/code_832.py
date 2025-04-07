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
    This function detects a strategy where a heterocycle is formed from precursors
    that already contain a functional linker that can be modified later.
    """
    # Track if we found the key features
    found_linker_installation = False
    found_heterocycle_formation = False
    found_linker_modification = False

    # Track the depth at which each feature was found
    linker_installation_depth = -1
    heterocycle_formation_depth = -1
    linker_modification_depth = -1

    # Track the molecules with installed linkers and heterocycles
    molecules_with_linkers = {}  # molecule SMILES -> set of linker FGs
    molecules_with_heterocycles = {}  # molecule SMILES -> set of heterocycles

    # List of heterocycles to check
    heterocycles = [
        "quinoline",
        "isoquinoline",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indole",
        "purine",
        "quinazoline",
        "pteridin",
        "imidazole",
        "oxazole",
        "thiazole",
        "furan",
        "thiophene",
        "pyrrole",
    ]

    # List of potential linker functional groups
    linker_fgs = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Ether",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Nitrile",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Phenol",
        "Aromatic halide",
    ]

    # List of linker installation reactions
    linker_installation_reactions = [
        "Williamson Ether Synthesis",
        "Williamson Ether Synthesis (intra to epoxy)",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "S-alkylation of thiols",
        "S-alkylation of thiols (ethyl)",
        "S-alkylation of thiols with alcohols",
        "S-alkylation of thiols with alcohols (ethyl)",
        "Esterification of Carboxylic Acids",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    ]

    # List of heterocycle formation reactions
    heterocycle_formation_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
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
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "Paal-Knorr pyrrole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
    ]

    # List of linker modification reactions
    linker_modification_reactions = [
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of alcohol to carboxylic acid",
        "Reduction of ester to primary alcohol",
        "Reduction of ketone to secondary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Alcohol to azide",
        "Nitrile to amide",
        "Alcohol to ether",
        "Ester saponification (methyl deprotection)",
        "Ester saponification (alkyl deprotection)",
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        "Reduction of primary amides to amines",
        "Reduction of secondary amides to amines",
        "Reduction of tertiary amides to amines",
        "Reduction of nitrile to amine",
        "Primary amine to fluoride",
        "Primary amine to chloride",
        "Primary amine to bromide",
        "Primary amine to iodide",
        "Alcohol to chloride_sulfonyl chloride",
        "Alcohol to chloride_SOCl2",
        "Alcohol to triflate conversion",
        "Appel reaction",
    ]

    # Track the synthesis path
    synthesis_path = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal found_linker_installation, found_heterocycle_formation, found_linker_modification
        nonlocal linker_installation_depth, heterocycle_formation_depth, linker_modification_depth
        nonlocal synthesis_path

        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for linker installation
                is_linker_installation = any(
                    checker.check_reaction(rxn, rsmi) for rxn in linker_installation_reactions
                )

                # If not a known linker installation reaction, check for FG changes that might indicate linker installation
                if not is_linker_installation:
                    # Check if any reactant has a potential linker attachment point
                    reactants_have_attachment_point = any(
                        checker.check_fg("Phenol", r)
                        or checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        or checker.check_fg("Carboxylic acid", r)
                        for r in reactants_smiles
                    )

                    # Check if product has a linker functional group that wasn't in the reactants
                    product_has_new_linker = False
                    for fg in linker_fgs:
                        if checker.check_fg(fg, product_smiles):
                            if not any(checker.check_fg(fg, r) for r in reactants_smiles):
                                product_has_new_linker = True
                                print(f"Found new linker FG: {fg} in product")
                                break

                    is_linker_installation = (
                        reactants_have_attachment_point and product_has_new_linker
                    )

                if is_linker_installation:
                    found_linker_installation = True
                    linker_installation_depth = depth

                    # Store which linker FGs are present in the product
                    linker_fgs_in_product = set()
                    for fg in linker_fgs:
                        if checker.check_fg(fg, product_smiles):
                            linker_fgs_in_product.add(fg)

                    molecules_with_linkers[product_smiles] = linker_fgs_in_product
                    print(f"Found linker installation at depth {depth}")
                    print(f"Linker FGs in product: {linker_fgs_in_product}")

                # Check for heterocycle formation
                is_heterocycle_formation = any(
                    checker.check_reaction(rxn, rsmi) for rxn in heterocycle_formation_reactions
                )

                # If not a known heterocycle formation reaction, check for ring formation
                if not is_heterocycle_formation:
                    # Check which heterocycles are in the product but not in any reactant
                    new_heterocycles = set()
                    for ring in heterocycles:
                        if checker.check_ring(ring, product_smiles):
                            if not any(checker.check_ring(ring, r) for r in reactants_smiles):
                                new_heterocycles.add(ring)
                                print(f"Found new heterocycle: {ring} in product")

                    is_heterocycle_formation = len(new_heterocycles) > 0

                if is_heterocycle_formation:
                    # Check if any reactant has a linker
                    reactants_with_linkers = []
                    for r in reactants_smiles:
                        if r in molecules_with_linkers:
                            reactants_with_linkers.append(r)

                    # Also check if any reactant has a linker functional group
                    reactants_with_linker_fgs = []
                    for r in reactants_smiles:
                        linker_fgs_in_reactant = set()
                        for fg in linker_fgs:
                            if checker.check_fg(fg, r):
                                linker_fgs_in_reactant.add(fg)
                        if linker_fgs_in_reactant:
                            reactants_with_linker_fgs.append((r, linker_fgs_in_reactant))

                    if reactants_with_linkers or reactants_with_linker_fgs:
                        # Check if the linker is preserved in the product
                        product_has_linker = False
                        for fg in linker_fgs:
                            if checker.check_fg(fg, product_smiles):
                                product_has_linker = True
                                break

                        if product_has_linker:
                            found_heterocycle_formation = True
                            heterocycle_formation_depth = depth

                            # Store which heterocycles are in the product
                            heterocycles_in_product = set()
                            for ring in heterocycles:
                                if checker.check_ring(ring, product_smiles):
                                    heterocycles_in_product.add(ring)

                            molecules_with_heterocycles[product_smiles] = heterocycles_in_product
                            print(
                                f"Found heterocycle formation with pre-installed linker at depth {depth}"
                            )
                            print(f"Heterocycles in product: {heterocycles_in_product}")

                # Check for linker modification after heterocycle formation
                if found_heterocycle_formation:
                    # Check if any reactant has a heterocycle
                    reactants_with_heterocycles = []
                    for r in reactants_smiles:
                        if r in molecules_with_heterocycles:
                            reactants_with_heterocycles.append(r)

                    # Also check if any reactant has a heterocycle directly
                    if not reactants_with_heterocycles:
                        for r in reactants_smiles:
                            heterocycles_in_reactant = set()
                            for ring in heterocycles:
                                if checker.check_ring(ring, r):
                                    heterocycles_in_reactant.add(ring)
                            if heterocycles_in_reactant:
                                reactants_with_heterocycles.append(r)

                    if reactants_with_heterocycles:
                        # Check if this is a known linker modification reaction
                        is_linker_modification = any(
                            checker.check_reaction(rxn, rsmi)
                            for rxn in linker_modification_reactions
                        )

                        # If not a known linker modification reaction, check for FG changes
                        if not is_linker_modification:
                            # Check if reactant has one FG and product has another
                            for fg in linker_fgs:
                                for other_fg in linker_fgs:
                                    if fg != other_fg:
                                        reactants_have_fg = any(
                                            checker.check_fg(fg, r) for r in reactants_smiles
                                        )
                                        product_has_other_fg = checker.check_fg(
                                            other_fg, product_smiles
                                        )

                                        if reactants_have_fg and product_has_other_fg:
                                            is_linker_modification = True
                                            print(f"Found FG change: {fg} -> {other_fg}")
                                            break
                                if is_linker_modification:
                                    break

                        # Verify the heterocycle is preserved in the product
                        product_has_heterocycle = False
                        for ring in heterocycles:
                            if checker.check_ring(ring, product_smiles):
                                product_has_heterocycle = True
                                break

                        if is_linker_modification and product_has_heterocycle:
                            found_linker_modification = True
                            linker_modification_depth = depth
                            print(f"Found linker modification at depth {depth}")

                            # Add this reaction to the synthesis path
                            synthesis_path.append(
                                {"type": "linker_modification", "depth": depth, "rsmi": rsmi}
                            )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # If this is a molecule node with a linker or heterocycle, track it
        elif node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for linker FGs
            linker_fgs_in_mol = set()
            for fg in linker_fgs:
                if checker.check_fg(fg, mol_smiles):
                    linker_fgs_in_mol.add(fg)

            if linker_fgs_in_mol:
                molecules_with_linkers[mol_smiles] = linker_fgs_in_mol

            # Check for heterocycles
            heterocycles_in_mol = set()
            for ring in heterocycles:
                if checker.check_ring(ring, mol_smiles):
                    heterocycles_in_mol.add(ring)

            if heterocycles_in_mol:
                molecules_with_heterocycles[mol_smiles] = heterocycles_in_mol

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Check if the key features were found in the correct order
    if found_linker_installation and found_heterocycle_formation:
        # In retrosynthetic traversal, higher depth means earlier in the synthesis
        if linker_installation_depth > heterocycle_formation_depth:
            print("Found strategy: linker installation followed by heterocycle formation")
            # If there was also a linker modification after heterocycle formation, that's a bonus
            if (
                found_linker_modification
                and linker_modification_depth < heterocycle_formation_depth
            ):
                print(
                    "Complete strategy found: linker installation → heterocycle formation → linker modification"
                )
            return True

    # If we didn't find the complete pattern in the correct order, check if we at least found
    # linker installation and heterocycle formation, even if the depths aren't in the expected order
    if found_linker_installation and found_heterocycle_formation:
        print("Found linker installation and heterocycle formation, but not in the expected order")
        return True

    return False
