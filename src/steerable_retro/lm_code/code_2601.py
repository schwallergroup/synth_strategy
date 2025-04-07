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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects a strategy involving sequential aromatic elaboration
    with nitro and methoxy groups, followed by carboxylic acid derivatization
    and late-stage benzylic functionalization.
    """
    # Track key transformations
    has_nitro_installation = False
    has_methoxy_installation = False
    has_carboxylic_derivatization = False
    has_benzylic_functionalization = False

    # Track reaction depths for determining sequence
    nitro_depth = -1
    methoxy_depth = -1
    carboxylic_depth = -1
    benzylic_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_installation, has_methoxy_installation, has_carboxylic_derivatization
        nonlocal has_benzylic_functionalization, nitro_depth, methoxy_depth, carboxylic_depth, benzylic_depth

        # Check if the final product already has the required functional groups
        if depth == 0 and node["type"] == "mol":
            product_smiles = node["smiles"]
            if not has_nitro_installation and checker.check_fg("Nitro group", product_smiles):
                has_nitro_installation = True
                nitro_depth = 999  # High depth to indicate it's present from the beginning
                print(f"Detected nitro group already present in final product")

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitro group installation
                if not has_nitro_installation:
                    # Check for nitration reactions
                    if (
                        checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                        or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                        or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                        or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                        or checker.check_reaction("Non-aromatic nitration with HNO3", rsmi)
                    ):
                        has_nitro_installation = True
                        nitro_depth = depth
                        print(f"Detected nitro group installation via reaction at depth {depth}")
                    else:
                        # Fallback to functional group check
                        product_has_nitro = checker.check_fg("Nitro group", product_smiles)
                        reactants_have_nitro = any(
                            checker.check_fg("Nitro group", r) for r in reactants_smiles
                        )

                        if product_has_nitro and not reactants_have_nitro:
                            has_nitro_installation = True
                            nitro_depth = depth
                            print(f"Detected nitro group installation at depth {depth}")

                # Check for methoxy group installation
                if not has_methoxy_installation:
                    # Check for methylation reactions
                    if (
                        checker.check_reaction("Methylation with MeI_primary", rsmi)
                        or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                        or checker.check_reaction("Methylation with MeI_aryl", rsmi)
                        or checker.check_reaction("Methylation with DMS", rsmi)
                        or checker.check_reaction("Methylation of OH with DMS", rsmi)
                        or checker.check_reaction("O-methylation", rsmi)
                        or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    ):
                        # Verify it's a methoxy on an aromatic ring
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            methoxy_pattern = Chem.MolFromSmarts("c-[O;X2][C;X4;H3]")
                            if product_mol.HasSubstructMatch(methoxy_pattern):
                                has_methoxy_installation = True
                                methoxy_depth = depth
                                print(
                                    f"Detected methoxy group installation via reaction at depth {depth}"
                                )
                    else:
                        # Fallback to functional group check
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

                        if product_mol and all(r for r in reactant_mols):
                            methoxy_pattern = Chem.MolFromSmarts("c-[O;X2][C;X4;H3]")
                            product_has_methoxy = product_mol.HasSubstructMatch(methoxy_pattern)
                            reactants_have_methoxy = any(
                                r.HasSubstructMatch(methoxy_pattern) for r in reactant_mols if r
                            )

                            if product_has_methoxy and not reactants_have_methoxy:
                                has_methoxy_installation = True
                                methoxy_depth = depth
                                print(f"Detected methoxy group installation at depth {depth}")

                # Check for carboxylic acid derivatization (including esterification)
                if not has_carboxylic_derivatization:
                    # Check for esterification reactions
                    if (
                        checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                        or checker.check_reaction("Transesterification", rsmi)
                        or checker.check_reaction("Acetic anhydride and alcohol to ester", rsmi)
                        or checker.check_reaction(
                            "O-alkylation of carboxylic acids with diazo compounds", rsmi
                        )
                        or checker.check_reaction("Carboxylic acid to amide conversion", rsmi)
                    ):
                        has_carboxylic_derivatization = True
                        carboxylic_depth = depth
                        print(f"Detected esterification/amidation at depth {depth}")
                    else:
                        # Check for acid to ester conversion by functional group change
                        product_has_ester = checker.check_fg("Ester", product_smiles)
                        product_has_amide = (
                            checker.check_fg("Primary amide", product_smiles)
                            or checker.check_fg("Secondary amide", product_smiles)
                            or checker.check_fg("Tertiary amide", product_smiles)
                        )
                        product_has_acyl = checker.check_fg("Acyl halide", product_smiles)
                        reactants_have_acid = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                        )

                        if (
                            product_has_ester or product_has_amide or product_has_acyl
                        ) and reactants_have_acid:
                            has_carboxylic_derivatization = True
                            carboxylic_depth = depth
                            print(f"Detected carboxylic acid derivatization at depth {depth}")

                # Check for benzylic functionalization
                if not has_benzylic_functionalization:
                    # Check for specific benzylic functionalization reactions
                    if any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Wohl-Ziegler bromination benzyl primary",
                            "Wohl-Ziegler bromination benzyl secondary",
                            "Wohl-Ziegler bromination benzyl tertiary",
                            "beta C(sp3) arylation",
                        ]
                    ):
                        has_benzylic_functionalization = True
                        benzylic_depth = depth
                        print(f"Detected benzylic functionalization via reaction at depth {depth}")
                    else:
                        # Convert to RDKit molecules for more detailed analysis
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

                        if product_mol and all(r for r in reactant_mols):
                            # Check for benzylic position in product
                            benzylic_pattern = Chem.MolFromSmarts("[c]-[C;!$(C=*);!$(C#*)]")

                            if product_mol.HasSubstructMatch(benzylic_pattern):
                                # Check if any reactant has a benzylic position that's being functionalized
                                for r_mol in reactant_mols:
                                    if r_mol.HasSubstructMatch(benzylic_pattern):
                                        # Compare benzylic positions in reactant and product
                                        product_matches = product_mol.GetSubstructMatches(
                                            benzylic_pattern
                                        )
                                        reactant_matches = r_mol.GetSubstructMatches(
                                            benzylic_pattern
                                        )

                                        # If there's a difference in the number of matches or their nature
                                        if len(product_matches) != len(reactant_matches):
                                            has_benzylic_functionalization = True
                                            benzylic_depth = depth
                                            print(
                                                f"Detected benzylic functionalization at depth {depth}"
                                            )
                                            break

                                        # Check for halogenation at benzylic position
                                        benzylic_halide = Chem.MolFromSmarts(
                                            "[c]-[C;!$(C=*);!$(C#*)]-[F,Cl,Br,I]"
                                        )
                                        if product_mol.HasSubstructMatch(
                                            benzylic_halide
                                        ) and not r_mol.HasSubstructMatch(benzylic_halide):
                                            has_benzylic_functionalization = True
                                            benzylic_depth = depth
                                            print(
                                                f"Detected benzylic halogenation at depth {depth}"
                                            )
                                            break
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present with correct sequence
    strategy_present = (
        has_nitro_installation
        and has_methoxy_installation
        and has_carboxylic_derivatization
        and has_benzylic_functionalization
        and
        # Benzylic functionalization should be late-stage (low depth)
        benzylic_depth <= 2
    )

    print(f"Strategy detected: {strategy_present}")
    print(f"Nitro depth: {nitro_depth}, Methoxy depth: {methoxy_depth}")
    print(f"Carboxylic depth: {carboxylic_depth}, Benzylic depth: {benzylic_depth}")

    return strategy_present
