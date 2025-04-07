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
    Detects the combined strategy of protection, cross-coupling, and SNAr in a single route.
    """
    has_protection = False
    has_cross_coupling = False
    has_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal has_protection, has_cross_coupling, has_snar

        indent = "  " * depth

        # Check if molecule has protected groups
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            print(f"{indent}Checking molecule: {mol_smiles}")

            # Check for Boc-protected amines
            if checker.check_fg("Carbamic ester", mol_smiles) and "C(C)(C)C" in mol_smiles:
                print(f"{indent}Found molecule with Boc group")
                has_protection = True

            # Check for silyl-protected alcohols
            if checker.check_fg("Silyl protective group", mol_smiles) or checker.check_fg(
                "TMS ether protective group", mol_smiles
            ):
                print(f"{indent}Found molecule with silyl protection")
                has_protection = True

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"{indent}Checking reaction: {rsmi}")

            # Check for protection reactions
            protection_reactions = [
                "Boc amine protection",
                "Boc amine protection explicit",
                "Boc amine protection with Boc anhydride",
                "Boc amine protection (ethyl Boc)",
                "Boc amine protection of secondary amine",
                "Boc amine protection of primary amine",
                "Alcohol protection with silyl ethers",
                "Protection of carboxylic acid",
            ]

            # Check for deprotection reactions (also indicates protection strategy)
            deprotection_reactions = [
                "Boc amine deprotection",
                "Alcohol deprotection from silyl ethers",
                "Alcohol deprotection from silyl ethers (double)",
                "Alcohol deprotection from silyl ethers (diol)",
                "Ester saponification (methyl deprotection)",
                "Ester saponification (alkyl deprotection)",
                "Hydroxyl benzyl deprotection",
                "Carboxyl benzyl deprotection",
                "Cleavage of methoxy ethers to alcohols",
                "Cleavage of alkoxy ethers to alcohols",
            ]

            # Check for protection reactions
            for rxn_name in protection_reactions:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"{indent}Found protection reaction: {rxn_name}")
                    has_protection = True
                    break

            # Check for deprotection reactions
            if not has_protection:
                for rxn_name in deprotection_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"{indent}Found deprotection reaction: {rxn_name}")
                        has_protection = True
                        break

            # Check for Boc protection specifically
            if not has_protection:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection
                if any("C(=O)OC(C)(C)C" in r for r in reactants) or "OC(=O)OC(C)(C)C" in "".join(
                    reactants
                ):
                    if "[NH" in "".join(reactants) and "NC(=O)OC(C)(C)C" in product:
                        print(f"{indent}Found Boc protection reaction by pattern matching")
                        has_protection = True

                # Check for amide formation (can be a form of protection)
                if "C(=O)" in product and "[NH" in "".join(reactants) and "NC(=O)" in product:
                    if not any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        ]
                    ):
                        print(f"{indent}Found potential amide protection by pattern matching")
                        has_protection = True

            # Check for cross-coupling reactions
            coupling_reactions = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with boronic esters",
                "Negishi coupling",
                "Stille reaction_vinyl",
                "Stille reaction_aryl",
                "Hiyama-Denmark Coupling",
                "Kumada cross-coupling",
                "Suzuki",
            ]

            for rxn_name in coupling_reactions:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"{indent}Found cross-coupling reaction: {rxn_name}")
                    has_cross_coupling = True
                    break

            # Check for Suzuki coupling specifically if not already found
            if not has_cross_coupling:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if boronic acid/ester and halide in reactants form a new C-C bond in product
                boronic_present = any("OB" in r for r in reactants)
                halide_present = any(("I" in r or "Br" in r or "Cl" in r) for r in reactants)

                if boronic_present and halide_present:
                    print(f"{indent}Found Suzuki coupling by pattern matching")
                    has_cross_coupling = True

            # Check for SNAr reactions
            snar_reactions = [
                "Ullmann-Goldberg Substitution amine",
                "Ullmann-Goldberg Substitution thiol",
                "Ullmann-Goldberg Substitution aryl alcohol",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "heteroaromatic_nuc_sub",
                "nucl_sub_aromatic_ortho_nitro",
                "nucl_sub_aromatic_para_nitro",
                "Buchwald-Hartwig",
                "N-arylation_heterocycles",
            ]

            for rxn_name in snar_reactions:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"{indent}Found SNAr reaction: {rxn_name}")
                    has_snar = True
                    break

            # Check for SNAr specifically if not already found
            if not has_snar:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if halide on aromatic ring is replaced by N, O, or S nucleophile
                halide_on_aromatic = any(
                    ("Cl" in r or "F" in r or "Br" in r or "I" in r) and ("c" in r or "n" in r)
                    for r in reactants
                )
                nucleophile_present = any(
                    ("[NH" in r or "[OH" in r or "[SH" in r or "NH3" in r) for r in reactants
                )

                if halide_on_aromatic and nucleophile_present:
                    print(f"{indent}Found SNAr by pattern matching")
                    has_snar = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Protection found: {has_protection}")
    print(f"Cross-coupling found: {has_cross_coupling}")
    print(f"SNAr found: {has_snar}")

    return has_protection and has_cross_coupling and has_snar
