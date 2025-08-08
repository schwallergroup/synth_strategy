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
    This function detects a strategy where a haloalkoxy linker is installed on an aromatic core,
    followed by scaffold-building reactions, and finally a late-stage amine substitution.
    """
    # Track if we found the key features
    found_linker_installation = False
    found_heterocycle_formation = False
    found_late_stage_amine_substitution = False

    # Track the depth at which each feature was found
    linker_installation_depth = -1
    heterocycle_formation_depth = -1
    amine_substitution_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_linker_installation, found_heterocycle_formation, found_late_stage_amine_substitution
        nonlocal linker_installation_depth, heterocycle_formation_depth, amine_substitution_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for linker installation (broader definition)
                # Linker installation can be via ether formation, ester formation, or amide formation
                if (
                    (
                        any(checker.check_fg("Phenol", r) for r in reactants_smiles)
                        and (
                            any(
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Aromatic halide", r)
                                or checker.check_fg("Triflate", r)
                                or checker.check_fg("Mesylate", r)
                                or checker.check_fg("Tosylate", r)
                                for r in reactants_smiles
                            )
                        )
                        and (
                            checker.check_reaction("Williamson Ether Synthesis", rsmi)
                            or checker.check_reaction("Mitsunobu_phenole", rsmi)
                            or checker.check_reaction("Williamson ether", rsmi)
                        )
                    )
                    or (
                        any(checker.check_fg("Phenol", r) for r in reactants_smiles)
                        and checker.check_fg("Ether", product_smiles)
                        and not all(checker.check_fg("Ether", r) for r in reactants_smiles)
                    )
                    or (
                        any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                        and any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            for r in reactants_smiles
                        )
                        and checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    )
                    or (
                        any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                        and any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants_smiles
                        )
                        and (
                            checker.check_reaction(
                                "Carboxylic acid with primary amine to amide", rsmi
                            )
                            or checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                            )
                        )
                    )
                ):
                    found_linker_installation = True
                    linker_installation_depth = depth
                    print(f"Found linker installation at depth {depth}, reaction: {rsmi}")

                # Check for heterocycle formation
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

                if all([product_mol] + reactant_mols):  # Ensure all molecules parsed correctly
                    product_ring_count = len(Chem.GetSSSR(product_mol))
                    reactant_ring_count = sum(len(Chem.GetSSSR(r)) for r in reactant_mols)

                    if product_ring_count > reactant_ring_count:
                        # Check if any of the new rings are heterocycles
                        heterocycle_rings = [
                            "pyridine",
                            "pyrimidine",
                            "pyrazine",
                            "pyridazine",
                            "triazine",
                            "quinoline",
                            "isoquinoline",
                            "quinazoline",
                            "quinoxaline",
                            "phthalazine",
                            "furan",
                            "thiophene",
                            "pyrrole",
                            "oxazole",
                            "thiazole",
                            "imidazole",
                            "triazole",
                            "tetrazole",
                            "piperidine",
                            "piperazine",
                            "morpholine",
                            "indole",
                            "benzimidazole",
                            "benzoxazole",
                            "benzothiazole",
                            "purine",
                        ]

                        for ring in heterocycle_rings:
                            if checker.check_ring(ring, product_smiles) and not any(
                                checker.check_ring(ring, r) for r in reactants_smiles
                            ):
                                found_heterocycle_formation = True
                                heterocycle_formation_depth = depth
                                print(f"Found heterocycle formation ({ring}) at depth {depth}")
                                break

                # Check for late-stage amine substitution (broader definition)
                amine_reactions = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "reductive amination",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "Goldberg coupling",
                    "Ullmann-Goldberg Substitution amine",
                    "aza-Michael addition primary",
                    "aza-Michael addition secondary",
                    "aza-Michael addition aromatic",
                    "Ring opening of epoxide with amine",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in amine_reactions):
                    found_late_stage_amine_substitution = True
                    amine_substitution_depth = depth
                    print(f"Found late-stage amine substitution at depth {depth}, reaction: {rsmi}")
                elif (
                    checker.check_fg("Secondary amine", product_smiles)
                    or checker.check_fg("Tertiary amine", product_smiles)
                ) and any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    or checker.check_fg("Aromatic halide", r)
                    or checker.check_fg("Aldehyde", r)
                    or checker.check_fg("Ketone", r)
                    or checker.check_fg("Triflate", r)
                    or checker.check_fg("Mesylate", r)
                    or checker.check_fg("Tosylate", r)
                    for r in reactants_smiles
                ):
                    found_late_stage_amine_substitution = True
                    amine_substitution_depth = depth
                    print(
                        f"Found late-stage amine substitution at depth {depth} based on functional groups"
                    )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all features were found
    print(f"Linker installation: {found_linker_installation} at depth {linker_installation_depth}")
    print(
        f"Heterocycle formation: {found_heterocycle_formation} at depth {heterocycle_formation_depth}"
    )
    print(
        f"Amine substitution: {found_late_stage_amine_substitution} at depth {amine_substitution_depth}"
    )

    # Check for the strategy with correct ordering
    # In retrosynthetic traversal, lower depth means later stage in synthesis
    if found_heterocycle_formation:  # At minimum, we need heterocycle formation
        # If we found all three features
        if found_linker_installation and found_late_stage_amine_substitution:
            # Linker installation should be at higher depth (earlier stage) than heterocycle formation
            # Amine substitution should be at lower depth (later stage) than heterocycle formation
            if (
                linker_installation_depth > heterocycle_formation_depth
                and heterocycle_formation_depth > amine_substitution_depth
            ):
                print(
                    "Found complete strategy: linker installation → heterocycle formation → late-stage amine substitution"
                )
                return True
        # If we only found heterocycle formation and amine substitution
        elif found_late_stage_amine_substitution:
            if heterocycle_formation_depth > amine_substitution_depth:
                print(
                    "Found partial strategy: heterocycle formation → late-stage amine substitution"
                )
                return True
        # If we only found heterocycle formation and linker installation
        elif found_linker_installation:
            if linker_installation_depth > heterocycle_formation_depth:
                print("Found partial strategy: linker installation → heterocycle formation")
                return True

    # If we found linker installation and amine substitution but no heterocycle formation
    # This is a special case that might still be valid
    if (
        found_linker_installation
        and found_late_stage_amine_substitution
        and not found_heterocycle_formation
    ):
        if linker_installation_depth > amine_substitution_depth:
            print("Found alternative strategy: linker installation → late-stage amine substitution")
            return True

    return False
