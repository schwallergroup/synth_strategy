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
    This function detects if the synthetic route involves modification of an amino acid
    side chain while preserving the amino acid backbone.
    """
    # Track amino acid backbone and modifications
    amino_acid_backbone_present = False
    side_chain_modifications = 0

    def dfs_traverse(node):
        nonlocal amino_acid_backbone_present, side_chain_modifications

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants_smiles = reactants_part.split(".")

                # Check for amino acid backbone in reactants and product
                reactant_has_aa = any(has_amino_acid_backbone(r) for r in reactants_smiles)
                product_has_aa = has_amino_acid_backbone(product_part)

                if reactant_has_aa and product_has_aa:
                    amino_acid_backbone_present = True
                    print(f"Found amino acid backbone in reaction: {rsmi}")

                    # Check for side chain modifications
                    if is_side_chain_modification(rsmi):
                        side_chain_modifications += 1
                        print(f"Found side chain modification: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    def has_amino_acid_backbone(smiles):
        """Check if molecule has an amino acid backbone (amine + carboxylic acid in correct arrangement)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for standard amino acid backbone
        # N-C-C(=O)-O or protected versions
        patterns = [
            # Standard amino acid: N-C-C(=O)-O
            "[NX3;!$(NC=O)]C([*])C(=O)[O,N]",
            # N-protected amino acid: C(=O)-N-C-C(=O)-O
            "[CX3](=O)[NX3]C([*])C(=O)[O,N]",
            # N-Boc protected: O=C(O[C])[NX3]C([*])C(=O)[O,N]
            "O=C(O[C])[NX3]C([*])C(=O)[O,N]",
            # N-Fmoc protected: O=C(OCC1c2ccccc2-c2ccccc21)[NX3]C([*])C(=O)[O,N]",
            "O=C(OCC1c2ccccc2-c2ccccc21)[NX3]C([*])C(=O)[O,N]",
            # C-protected (ester): N-C-C(=O)-OC
            "[NX3;!$(NC=O)]C([*])C(=O)O[C]",
            # C-protected (amide): N-C-C(=O)-N
            "[NX3;!$(NC=O)]C([*])C(=O)[NX3]",
        ]

        for pattern in patterns:
            patt = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(patt):
                return True

        return False

    def is_side_chain_modification(rsmi):
        """Check if the reaction modifies the side chain while preserving the backbone"""
        reactants_part = rsmi.split(">")[0]
        product_part = rsmi.split(">")[-1]

        # Check for common side chain modification reactions
        common_modifications = [
            "Oxidation of aldehydes to carboxylic acids",
            "Oxidation of alcohols to carboxylic acids",
            "Oxidation of primary alcohols",
            "Reduction of carboxylic acid to primary alcohol",
            "Reduction of ester to primary alcohol",
            "Esterification of Carboxylic Acids",
            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
            "Protection of carboxylic acid",
            "Deprotection of carboxylic acid",
            "Boc amine protection",
            "Boc amine deprotection",
            "Alcohol protection with silyl ethers",
            "Alcohol deprotection from silyl ethers",
            "Williamson Ether Synthesis",
            "Alkylation of amines",
            "Reductive amination with aldehyde",
            "Reductive amination with ketone",
            "Acylation of primary amines",
            "Acylation of secondary amines",
            "Sulfonamide synthesis (Schotten-Baumann) primary amine",
            "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        ]

        for rxn_type in common_modifications:
            if checker.check_reaction(rxn_type, rsmi):
                # Verify the amino acid backbone is preserved
                if backbone_preserved(reactants_part, product_part):
                    return True

        # Check for functional group transformations on side chains
        reactants = reactants_part.split(".")
        reactant_with_aa = None

        # Find the reactant with amino acid backbone
        for r in reactants:
            if has_amino_acid_backbone(r):
                reactant_with_aa = r
                break

        if not reactant_with_aa:
            return False

        # Check for changes in functional groups that might be on side chains
        fg_types = [
            "Primary alcohol",
            "Secondary alcohol",
            "Tertiary alcohol",
            "Carboxylic acid",
            "Ester",
            "Primary amide",
            "Secondary amide",
            "Primary amine",
            "Secondary amine",
            "Tertiary amine",
            "Thiol",
            "Thioester",
            "Nitrile",
            "Aldehyde",
            "Ketone",
            "Phenol",
            "Aromatic halide",
            "Primary halide",
            "Secondary halide",
            "Tertiary halide",
            "Azide",
            "Alkyne",
            "Aromatic alcohol",
            "Sulfonamide",
            "Sulfone",
            "Sulfoxide",
        ]

        for fg in fg_types:
            reactant_has_fg = checker.check_fg(fg, reactant_with_aa)
            product_has_fg = checker.check_fg(fg, product_part)

            if reactant_has_fg != product_has_fg:
                # Verify the amino acid backbone is preserved
                if backbone_preserved(reactants_part, product_part):
                    return True

        return False

    def backbone_preserved(reactants_smiles, product_smiles):
        """Check if the amino acid backbone is preserved between reactants and product"""
        # Check for amino acid patterns in both reactants and product
        patterns = [
            # Standard amino acid: N-C-C(=O)-O
            "[NX3;!$(NC=O)]C([*])C(=O)[O,N]",
            # N-protected amino acid: C(=O)-N-C-C(=O)-O
            "[CX3](=O)[NX3]C([*])C(=O)[O,N]",
            # N-Boc protected: O=C(O[C])[NX3]C([*])C(=O)[O,N]
            "O=C(O[C])[NX3]C([*])C(=O)[O,N]",
            # N-Fmoc protected: O=C(OCC1c2ccccc2-c2ccccc21)[NX3]C([*])C(=O)[O,N]",
            "O=C(OCC1c2ccccc2-c2ccccc21)[NX3]C([*])C(=O)[O,N]",
            # C-protected (ester): N-C-C(=O)-OC
            "[NX3;!$(NC=O)]C([*])C(=O)O[C]",
            # C-protected (amide): N-C-C(=O)-N
            "[NX3;!$(NC=O)]C([*])C(=O)[NX3]",
        ]

        reactants_list = reactants_smiles.split(".")
        reactant_with_aa = None

        # Find reactant with amino acid backbone
        for r in reactants_list:
            mol = Chem.MolFromSmiles(r)
            if not mol:
                continue

            for pattern in patterns:
                patt = Chem.MolFromSmarts(pattern)
                if mol.HasSubstructMatch(patt):
                    reactant_with_aa = r
                    break
            if reactant_with_aa:
                break

        if not reactant_with_aa:
            return False

        # Check if product has amino acid backbone
        product_mol = Chem.MolFromSmiles(product_smiles)
        if not product_mol:
            return False

        for pattern in patterns:
            patt = Chem.MolFromSmarts(pattern)
            if product_mol.HasSubstructMatch(patt):
                return True

        return False

    # Start traversal
    dfs_traverse(route)

    # Check if amino acid side chain modification strategy was used
    result = amino_acid_backbone_present and side_chain_modifications >= 1
    print(f"Amino acid side chain modification strategy detected: {result}")
    print(f"  Amino acid backbone present: {amino_acid_backbone_present}")
    print(f"  Number of side chain modifications: {side_chain_modifications}")

    return result
