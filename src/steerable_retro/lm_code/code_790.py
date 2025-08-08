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
    This function detects a strategy involving late-stage N-methylation of a nitrogen heterocycle.
    """
    found_n_methylation = False

    # List of nitrogen heterocycles to check
    nitrogen_heterocycles = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "acridine",
        "benzimidazole",
        "indazole",
        "benzotriazole",
    ]

    # List of N-methylation reaction types
    methylation_reactions = [
        "N-methylation",
        "Eschweiler-Clarke Primary Amine Methylation",
        "Eschweiler-Clarke Secondary Amine Methylation",
        "Reductive methylation of primary amine with formaldehyde",
        "Methylation with MeI_primary",
        "Methylation with MeI_secondary",
        "Methylation with MeI_tertiary",
        "Methylation with DMS",
        "DMS Amine methylation",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_n_methylation

        if node["type"] == "reaction" and depth <= 3:  # Late stage (final 4 steps)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an N-methylation reaction
                is_methylation_reaction = False
                for rxn_type in methylation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected {rxn_type} reaction")
                        is_methylation_reaction = True
                        break

                # Additional check for N-methylation pattern in SMILES
                if not is_methylation_reaction:
                    # Check if product has a methyl group attached to nitrogen
                    # and if any reactant has an NH group (indicating potential methylation)
                    if ("[CH3][N" in product or "[CH3:1][N" in product) and any(
                        "[NH" in reactant or "[NH:" in reactant
                        for reactant in reactants
                        if reactant.strip()
                    ):
                        print("Detected N-methylation pattern in SMILES")
                        is_methylation_reaction = True

                    # Check for formaldehyde (common methylating agent) in reactants
                    if any(
                        "O=[CH2]" in reactant or "O=[CH2:" in reactant
                        for reactant in reactants
                        if reactant.strip()
                    ):
                        print("Detected formaldehyde (potential methylating agent)")
                        is_methylation_reaction = True

                if not is_methylation_reaction:
                    print("Not an N-methylation reaction")
                    for child in node.get("children", []):
                        dfs_traverse(child, depth + 1)
                    return

                # Check if product contains a nitrogen heterocycle
                product_has_n_heterocycle = False
                product_heterocycle = None
                for ring in nitrogen_heterocycles:
                    if checker.check_ring(ring, product):
                        print(f"Product contains {ring} ring")
                        product_has_n_heterocycle = True
                        product_heterocycle = ring
                        break

                if not product_has_n_heterocycle:
                    print("Product does not contain a nitrogen heterocycle")
                    for child in node.get("children", []):
                        dfs_traverse(child, depth + 1)
                    return

                # Check if reactants also contain the same nitrogen heterocycle
                reactant_has_same_heterocycle = False
                for reactant in reactants:
                    if len(reactant.strip()) > 0:  # Skip empty reactants
                        if checker.check_ring(product_heterocycle, reactant):
                            reactant_has_same_heterocycle = True
                            print(f"Reactant contains the same {product_heterocycle} ring")
                            break

                if not reactant_has_same_heterocycle:
                    print("No reactant contains the same heterocycle as the product")
                    for child in node.get("children", []):
                        dfs_traverse(child, depth + 1)
                    return

                # Check if reactants had primary or secondary amines that could be methylated
                reactant_has_methylatable_n = False
                for reactant in reactants:
                    if len(reactant.strip()) > 0:  # Skip empty reactants
                        if checker.check_ring(product_heterocycle, reactant):
                            if (
                                checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Primary amine", reactant)
                                or "[NH]" in reactant
                                or "[NH:" in reactant
                            ):
                                print(f"Reactant has methylatable nitrogen: {reactant}")
                                reactant_has_methylatable_n = True
                                break

                # Check if product has tertiary amine (result of methylation)
                has_tertiary_amine = checker.check_fg("Tertiary amine", product)
                print(f"Product has tertiary amine: {has_tertiary_amine}")

                # Additional check for methyl group on nitrogen in product
                has_n_methyl = "[CH3][N" in product or "[CH3:1][N" in product
                print(f"Product has N-methyl group: {has_n_methyl}")

                if reactant_has_methylatable_n and (has_tertiary_amine or has_n_methyl):
                    print(
                        f"Found late-stage N-methylation of a nitrogen heterocycle at depth {depth}"
                    )
                    found_n_methylation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_n_methylation
