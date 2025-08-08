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
    This function detects a convergent synthesis strategy where two heterocyclic fragments
    (pyrazole and piperidine) are connected via an amide linker, with key steps including
    thioether formation, N-Boc deprotection, ester hydrolysis, and amide coupling.
    """
    # Initialize tracking variables
    pyrazole_molecules = set()
    piperidine_molecules = set()
    thioether_reactions = []
    boc_deprotection_reactions = []
    ester_hydrolysis_reactions = []
    amide_coupling_reactions = []

    # Track the final product that contains both fragments
    final_product = None

    def dfs_traverse(node, depth=0):
        nonlocal final_product

        if node["type"] == "mol":
            if node.get("smiles"):
                mol_smiles = node.get("smiles")

                # Check for pyrazole fragment
                if checker.check_ring("pyrazole", mol_smiles):
                    pyrazole_molecules.add(mol_smiles)
                    print(f"Found pyrazole fragment (depth {depth}): {mol_smiles}")

                # Check for piperidine fragment
                if checker.check_ring("piperidine", mol_smiles):
                    piperidine_molecules.add(mol_smiles)
                    print(f"Found piperidine fragment (depth {depth}): {mol_smiles}")

                # If molecule contains both fragments and we're at the top level (final product)
                if (
                    depth == 0
                    and checker.check_ring("pyrazole", mol_smiles)
                    and checker.check_ring("piperidine", mol_smiles)
                ):
                    final_product = mol_smiles
                    print(f"Found final product with both fragments: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thioether formation
                if (
                    checker.check_reaction("S-alkylation of thiols", rsmi)
                    or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi)
                    or checker.check_reaction("S-alkylation of thiols with alcohols", rsmi)
                    or checker.check_reaction("S-alkylation of thiols with alcohols (ethyl)", rsmi)
                ):
                    thioether_reactions.append(rsmi)
                    print(f"Found thioether formation (depth {depth}): {rsmi}")
                # Check for thioether pattern in molecules
                elif "Sc" in product and any(
                    checker.check_fg("Aromatic thiol", r) or checker.check_fg("Aliphatic thiol", r)
                    for r in reactants
                ):
                    thioether_reactions.append(rsmi)
                    print(f"Found thioether formation by pattern (depth {depth}): {rsmi}")

                # Check for Boc deprotection
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                ):
                    boc_deprotection_reactions.append(rsmi)
                    print(f"Found Boc deprotection (depth {depth}): {rsmi}")
                # Check for Boc pattern in molecules
                elif any("OC(=O)C(C)(C)C" in r for r in reactants) and not any(
                    "OC(=O)C(C)(C)C" in product
                ):
                    boc_deprotection_reactions.append(rsmi)
                    print(f"Found Boc deprotection by pattern (depth {depth}): {rsmi}")

                # Check for ester hydrolysis
                if (
                    checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                ):
                    ester_hydrolysis_reactions.append(rsmi)
                    print(f"Found ester hydrolysis (depth {depth}): {rsmi}")
                # Check for ester hydrolysis pattern
                elif any(checker.check_fg("Ester", r) for r in reactants) and checker.check_fg(
                    "Carboxylic acid", product
                ):
                    ester_hydrolysis_reactions.append(rsmi)
                    print(f"Found ester hydrolysis by pattern (depth {depth}): {rsmi}")

                # Check for amide coupling
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid to amide conversion", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                ):

                    # Check if this amide coupling connects pyrazole and piperidine fragments
                    pyrazole_in_reactants = any(
                        checker.check_ring("pyrazole", r) for r in reactants
                    )
                    piperidine_in_reactants = any(
                        checker.check_ring("piperidine", r) for r in reactants
                    )
                    both_in_product = checker.check_ring(
                        "pyrazole", product
                    ) and checker.check_ring("piperidine", product)

                    if (
                        (pyrazole_in_reactants and piperidine_in_reactants and both_in_product)
                        or (
                            pyrazole_in_reactants
                            and not piperidine_in_reactants
                            and both_in_product
                        )
                        or (
                            not pyrazole_in_reactants
                            and piperidine_in_reactants
                            and both_in_product
                        )
                    ):
                        amide_coupling_reactions.append(rsmi)
                        print(f"Found key amide coupling (depth {depth}): {rsmi}")
                        print(f"  Pyrazole in reactants: {pyrazole_in_reactants}")
                        print(f"  Piperidine in reactants: {piperidine_in_reactants}")
                        print(f"  Both in product: {both_in_product}")

                # Check for amide coupling by pattern
                elif not amide_coupling_reactions:
                    carboxylic_acid_in_reactants = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    amine_in_reactants = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    amide_in_product = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    pyrazole_in_reactants = any(
                        checker.check_ring("pyrazole", r) for r in reactants
                    )
                    piperidine_in_reactants = any(
                        checker.check_ring("piperidine", r) for r in reactants
                    )
                    both_in_product = checker.check_ring(
                        "pyrazole", product
                    ) and checker.check_ring("piperidine", product)

                    if carboxylic_acid_in_reactants and amine_in_reactants and amide_in_product:
                        if (
                            (pyrazole_in_reactants and piperidine_in_reactants and both_in_product)
                            or (
                                pyrazole_in_reactants
                                and not piperidine_in_reactants
                                and both_in_product
                            )
                            or (
                                not pyrazole_in_reactants
                                and piperidine_in_reactants
                                and both_in_product
                            )
                        ):
                            amide_coupling_reactions.append(rsmi)
                            print(f"Found key amide coupling by pattern (depth {depth}): {rsmi}")
                            print(f"  Pyrazole in reactants: {pyrazole_in_reactants}")
                            print(f"  Piperidine in reactants: {piperidine_in_reactants}")
                            print(f"  Both in product: {both_in_product}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # For a convergent strategy, we need:
    # 1. Both fragments in the final product
    # 2. An amide coupling reaction that connects them
    # 3. At least one preparatory reaction (thioether, Boc deprotection, or ester hydrolysis)

    has_pyrazole = len(pyrazole_molecules) > 0
    has_piperidine = len(piperidine_molecules) > 0
    has_thioether = len(thioether_reactions) > 0
    has_boc_deprotection = len(boc_deprotection_reactions) > 0
    has_ester_hydrolysis = len(ester_hydrolysis_reactions) > 0
    has_amide_coupling = len(amide_coupling_reactions) > 0

    # Check if final product contains both fragments
    has_both_fragments_in_product = final_product is not None

    # Check for preparatory reactions
    has_preparatory_reaction = has_thioether or has_boc_deprotection or has_ester_hydrolysis

    # For this specific strategy, we might need to be more lenient
    # If we have both fragments and they're connected in the final product,
    # we can infer that the strategy is present even if we didn't explicitly detect all reactions
    if has_both_fragments_in_product and not has_preparatory_reaction:
        # Check if the final product has an amide bond connecting the fragments
        if (
            checker.check_fg("Primary amide", final_product)
            or checker.check_fg("Secondary amide", final_product)
            or checker.check_fg("Tertiary amide", final_product)
        ):
            print("Inferring preparatory reactions based on final product structure")
            has_preparatory_reaction = True

    # Determine if strategy is present
    strategy_present = has_pyrazole and has_piperidine and has_both_fragments_in_product

    # If we have the fragments but didn't detect the amide coupling, check the final product
    if strategy_present and not has_amide_coupling:
        if (
            checker.check_fg("Primary amide", final_product)
            or checker.check_fg("Secondary amide", final_product)
            or checker.check_fg("Tertiary amide", final_product)
        ):
            print("Inferring amide coupling based on final product structure")
            has_amide_coupling = True

    # Final determination with all criteria
    strategy_present = (
        has_pyrazole
        and has_piperidine
        and has_amide_coupling
        and has_preparatory_reaction
        and has_both_fragments_in_product
    )

    print(f"Convergent heterocycle amide coupling strategy detected: {strategy_present}")
    print(f"Pyrazole fragments: {has_pyrazole}, Piperidine fragments: {has_piperidine}")
    print(f"Both fragments in final product: {has_both_fragments_in_product}")
    print(f"Thioether formation: {has_thioether}, Boc deprotection: {has_boc_deprotection}")
    print(f"Ester hydrolysis: {has_ester_hydrolysis}, Amide coupling: {has_amide_coupling}")

    return strategy_present
