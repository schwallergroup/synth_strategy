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
    Detects if the synthesis includes protection with a benzyl group.
    """
    has_benzyl_protection = False

    def is_benzyl_halide(smiles):
        """Check if a molecule is a benzyl halide"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check if it's a halide
        if not (
            checker.check_fg("Primary halide", smiles)
            or checker.check_fg("Secondary halide", smiles)
        ):
            return False

        # Check for benzyl pattern
        benzyl_pattern = Chem.MolFromSmarts("c1ccccc1C[F,Cl,Br,I]")
        return mol.HasSubstructMatch(benzyl_pattern)

    def is_benzyl_alcohol(smiles):
        """Check if a molecule is benzyl alcohol or derivative"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        benzyl_alcohol_pattern = Chem.MolFromSmarts("c1ccccc1CO")
        return mol.HasSubstructMatch(benzyl_alcohol_pattern)

    def has_benzyl_ether(smiles):
        """Check if a molecule contains a benzyl ether group"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        benzyl_ether_pattern = Chem.MolFromSmarts("OCc1ccccc1")
        return mol.HasSubstructMatch(benzyl_ether_pattern)

    def has_benzyl_amine(smiles):
        """Check if a molecule contains a benzyl amine group"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        benzyl_amine_pattern = Chem.MolFromSmarts("NCc1ccccc1")
        return mol.HasSubstructMatch(benzyl_amine_pattern)

    def has_benzyl_ester(smiles):
        """Check if a molecule contains a benzyl ester group"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        benzyl_ester_pattern = Chem.MolFromSmarts("C(=O)OCc1ccccc1")
        return mol.HasSubstructMatch(benzyl_ester_pattern)

    def has_benzyl_thioether(smiles):
        """Check if a molecule contains a benzyl thioether group"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        benzyl_thioether_pattern = Chem.MolFromSmarts("SCc1ccccc1")
        return mol.HasSubstructMatch(benzyl_thioether_pattern)

    def dfs_traverse(node):
        nonlocal has_benzyl_protection

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for benzyl protection reactions
            # 1. Check for Williamson ether synthesis (O-benzylation)
            if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                print("Detected Williamson Ether Synthesis")
                for reactant in reactants:
                    if (
                        checker.check_fg("Phenol", reactant)
                        or checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                    ):
                        # Check if benzyl halide is present in reactants
                        for r in reactants:
                            if is_benzyl_halide(r):
                                print(
                                    f"Detected benzyl protection via Williamson ether synthesis: {rsmi}"
                                )
                                has_benzyl_protection = True

            # 2. Check for N-alkylation with benzyl halides
            if (
                checker.check_reaction("N-alkylation of primary amines with alkyl halides", rsmi)
                or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )
                or checker.check_reaction("Alkylation of amines", rsmi)
            ):
                print("Detected N-alkylation reaction")
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        # Check if benzyl halide is present in reactants
                        for r in reactants:
                            if is_benzyl_halide(r):
                                print(f"Detected benzyl protection via N-alkylation: {rsmi}")
                                has_benzyl_protection = True

            # 3. Check for esterification with benzyl alcohol
            if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                print("Detected Esterification reaction")
                for reactant in reactants:
                    if checker.check_fg("Carboxylic acid", reactant):
                        # Check if benzyl alcohol is present in reactants
                        for r in reactants:
                            if is_benzyl_alcohol(r):
                                print(f"Detected benzyl protection via esterification: {rsmi}")
                                has_benzyl_protection = True

            # 4. Check for S-alkylation with benzyl halides
            if (
                checker.check_reaction("S-alkylation of thiols", rsmi)
                or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi)
                or checker.check_reaction("S-alkylation of thiols with alcohols", rsmi)
                or checker.check_reaction("S-alkylation of thiols with alcohols (ethyl)", rsmi)
            ):
                print("Detected S-alkylation reaction")
                for reactant in reactants:
                    if checker.check_fg("Aromatic thiol", reactant) or checker.check_fg(
                        "Aliphatic thiol", reactant
                    ):
                        # Check if benzyl halide is present in reactants
                        for r in reactants:
                            if is_benzyl_halide(r) or is_benzyl_alcohol(r):
                                print(f"Detected benzyl protection via S-alkylation: {rsmi}")
                                has_benzyl_protection = True

            # 5. Check for benzyl deprotection reactions
            if (
                checker.check_reaction("Hydroxyl benzyl deprotection", rsmi)
                or checker.check_reaction("Carboxyl benzyl deprotection", rsmi)
                or checker.check_reaction("Hydrogenolysis of amides/imides/carbamates", rsmi)
                or checker.check_reaction("Hydrogenolysis of tertiary amines", rsmi)
            ):
                print(f"Detected benzyl deprotection reaction: {rsmi}")
                has_benzyl_protection = True

            # 6. Check for Mitsunobu reactions
            if (
                checker.check_reaction("Mitsunobu aryl ether", rsmi)
                or checker.check_reaction("Mitsunobu esterification", rsmi)
                or checker.check_reaction("Mitsunobu aryl ether (intramolecular)", rsmi)
                or checker.check_reaction("{Mitsunobu_imide}", rsmi)
                or checker.check_reaction("{Mitsunobu_phenole}", rsmi)
                or checker.check_reaction("{Mitsunobu_sulfonamide}", rsmi)
            ):
                print("Detected Mitsunobu reaction")
                for reactant in reactants:
                    if is_benzyl_alcohol(reactant):
                        print(f"Detected benzyl protection via Mitsunobu: {rsmi}")
                        has_benzyl_protection = True

            # 7. Check for direct benzyl protection patterns in product vs reactants
            for reactant in reactants:
                # Check for OH in reactant and O-Bn in product
                if (
                    checker.check_fg("Phenol", reactant)
                    or checker.check_fg("Primary alcohol", reactant)
                    or checker.check_fg("Secondary alcohol", reactant)
                    or checker.check_fg("Tertiary alcohol", reactant)
                ):
                    if has_benzyl_ether(product) and not has_benzyl_ether(reactant):
                        print(f"Detected benzyl protection of alcohol: {rsmi}")
                        has_benzyl_protection = True

                # Check for NH in reactant and N-Bn in product
                if (
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    or checker.check_fg("Aniline", reactant)
                ):
                    if has_benzyl_amine(product) and not has_benzyl_amine(reactant):
                        print(f"Detected benzyl protection of amine: {rsmi}")
                        has_benzyl_protection = True

                # Check for COOH in reactant and COO-Bn in product
                if checker.check_fg("Carboxylic acid", reactant):
                    if has_benzyl_ester(product) and not has_benzyl_ester(reactant):
                        print(f"Detected benzyl protection of carboxylic acid: {rsmi}")
                        has_benzyl_protection = True

                # Check for SH in reactant and S-Bn in product
                if checker.check_fg("Aromatic thiol", reactant) or checker.check_fg(
                    "Aliphatic thiol", reactant
                ):
                    if has_benzyl_thioether(product) and not has_benzyl_thioether(reactant):
                        print(f"Detected benzyl protection of thiol: {rsmi}")
                        has_benzyl_protection = True

            # 8. Check for benzyl group in product that wasn't in reactants (general case)
            benzyl_in_reactants = any(
                (
                    Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1C"))
                    if Chem.MolFromSmiles(r)
                    else False
                )
                for r in reactants
            )
            benzyl_in_product = (
                Chem.MolFromSmiles(product).HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1C"))
                if Chem.MolFromSmiles(product)
                else False
            )

            if benzyl_in_product and not benzyl_in_reactants:
                print(f"Detected addition of benzyl group: {rsmi}")
                has_benzyl_protection = True

        # Check if the molecule itself contains a protected benzyl group
        elif node["type"] == "mol" and node.get("smiles"):
            smiles = node["smiles"]
            if (
                has_benzyl_ether(smiles)
                or has_benzyl_amine(smiles)
                or has_benzyl_ester(smiles)
                or has_benzyl_thioether(smiles)
            ):
                print(f"Detected molecule with benzyl protection: {smiles}")
                has_benzyl_protection = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: has_benzyl_protection = {has_benzyl_protection}")
    return has_benzyl_protection
