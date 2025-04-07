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
    Detects a synthetic strategy involving late-stage introduction of
    cyclic amine fragments to a heterocyclic core.
    """
    # Track fragment introductions
    fragment_introductions = []
    depths = []

    # List of cyclic amines to check
    cyclic_amines = [
        "morpholine",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "thiomorpholine",
    ]

    # List of heterocyclic rings to check in the core
    heterocyclic_rings = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "indole",
        "benzimidazole",
        "quinoline",
        "isoquinoline",
        "triazole",
        "tetrazole",
        "benzoxazole",
        "benzothiazole",
        "purine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                # Extract depth safely
                depth = 0
                if "metadata" in node and "ID" in node["metadata"]:
                    id_str = node["metadata"]["ID"]
                    depth_match = re.search(r"Depth: (\d+)", id_str)
                    if depth_match:
                        depth = int(depth_match.group(1))

                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = reactants_str.split(".")

                # Check if product contains a heterocyclic core
                has_heterocyclic_core = False
                core_found = None
                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, product_str):
                        has_heterocyclic_core = True
                        core_found = ring
                        print(f"Found heterocyclic core: {ring} in product")
                        break

                if has_heterocyclic_core:
                    # Check for cyclic amine fragments in reactants
                    for reactant in reactants:
                        for amine in cyclic_amines:
                            if checker.check_ring(amine, reactant):
                                print(f"Found {amine} in reactant: {reactant}")

                                # Verify this is a fragment introduction by checking:
                                # 1. The amine is not in other reactants
                                # 2. The amine is in the product
                                other_reactants = [r for r in reactants if r != reactant]
                                other_reactants_have_amine = any(
                                    checker.check_ring(amine, r) for r in other_reactants
                                )

                                if not other_reactants_have_amine and checker.check_ring(
                                    amine, product_str
                                ):
                                    # Check if this is a nucleophilic substitution or coupling reaction
                                    is_coupling = (
                                        checker.check_reaction(
                                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                                            rsmi,
                                        )
                                        or checker.check_reaction(
                                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                            rsmi,
                                        )
                                        or checker.check_reaction(
                                            "Williamson Ether Synthesis", rsmi
                                        )
                                        or checker.check_reaction(
                                            "N-alkylation of primary amines with alkyl halides",
                                            rsmi,
                                        )
                                        or checker.check_reaction(
                                            "N-alkylation of secondary amines with alkyl halides",
                                            rsmi,
                                        )
                                        or checker.check_reaction("Buchwald-Hartwig", rsmi)
                                        or checker.check_reaction("N-arylation_heterocycles", rsmi)
                                        or checker.check_reaction(
                                            "Ullmann-Goldberg Substitution amine", rsmi
                                        )
                                        or checker.check_reaction("Goldberg coupling", rsmi)
                                        or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                                        or checker.check_reaction(
                                            "nucl_sub_aromatic_ortho_nitro", rsmi
                                        )
                                        or checker.check_reaction(
                                            "nucl_sub_aromatic_para_nitro", rsmi
                                        )
                                        or
                                        # Additional nucleophilic substitution reactions
                                        checker.check_reaction(
                                            "Ullmann-Goldberg Substitution aryl alcohol", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Goldberg coupling aryl amine-aryl chloride", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Goldberg coupling aryl amide-aryl chloride", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Displacement of ethoxy group by primary amine", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Displacement of ethoxy group by secondary amine", rsmi
                                        )
                                    )

                                    # Also check for general nucleophilic substitution patterns
                                    if not is_coupling:
                                        # Check if reactant has a halide (common leaving group)
                                        has_leaving_group = False
                                        for leaving_group in [
                                            "Primary halide",
                                            "Secondary halide",
                                            "Tertiary halide",
                                            "Aromatic halide",
                                        ]:
                                            for r in other_reactants:
                                                if checker.check_fg(leaving_group, r):
                                                    has_leaving_group = True
                                                    print(
                                                        f"Found leaving group: {leaving_group} in other reactant"
                                                    )
                                                    break
                                            if has_leaving_group:
                                                break

                                        # If there's a leaving group and the amine is nucleophilic, it's likely a substitution
                                        if has_leaving_group:
                                            is_coupling = True
                                            print(
                                                f"Inferred nucleophilic substitution based on leaving group"
                                            )

                                    if is_coupling:
                                        fragment_introductions.append((amine, depth))
                                        print(
                                            f"Detected {amine} introduction at depth {depth} with core {core_found}"
                                        )

                depths.append(depth)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Calculate max depth to determine what's "late stage"
    max_depth = max(depths) if depths else 0

    # Consider first third of synthesis as "late stage" (lower depth values)
    late_stage_threshold = max_depth // 3

    # Count late stage fragment introductions
    late_stage_introductions = sum(
        1 for _, depth in fragment_introductions if depth <= late_stage_threshold
    )

    # Determine if the strategy is present (at least 1 fragment introduced in late stage)
    strategy_present = late_stage_introductions >= 1

    print(f"Late-stage fragment introduction strategy detected: {strategy_present}")
    print(f"- Total fragment introductions: {len(fragment_introductions)}")
    print(f"- Late-stage introductions: {late_stage_introductions}")
    print(f"- Max depth: {max_depth}, Late stage threshold: {late_stage_threshold}")

    return strategy_present
