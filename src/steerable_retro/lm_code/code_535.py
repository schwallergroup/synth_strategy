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
    Detects if the route uses a late-stage fragment coupling strategy,
    specifically looking for multiple fragments combined in the final step.
    """
    late_coupling = False
    final_product_smiles = route["smiles"]
    print(f"Final product SMILES: {final_product_smiles}")

    # Important ring structures that should be considered substantial fragments
    important_rings = [
        "benzene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "indole",
        "quinoline",
        "isoquinoline",
        "naphthalene",
        "thiophene",
        "furan",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "piperidine",
        "piperazine",
        "morpholine",
        "tetrahydrofuran",
        "cyclohexane",
    ]

    # Important functional groups that indicate substantial fragments
    important_fgs = [
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Nitrile",
        "Sulfonamide",
        "Nitro group",
        "Boronic acid",
        "Boronic ester",
        "Aromatic halide",
    ]

    # Expanded list of coupling reactions
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with sulfonic esters",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_aryl OTf",
        "Stille reaction_vinyl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Stille reaction_other",
        "Stille reaction_other OTf",
        "Negishi coupling",
        "Heck terminal vinyl",
        "Heck_terminal_vinyl",
        "Heck_non-terminal_vinyl",
        "Oxidative Heck reaction",
        "Oxidative Heck reaction with vinyl ester",
        "Heck reaction with vinyl ester and amine",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl OTf",
        "Sonogashira acetylene_aryl OTf",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl halide",
        "Sonogashira alkyne_alkenyl OTf",
        "Sonogashira acetylene_alkenyl OTf",
        "Sonogashira alkyne_acyl halide",
        "Sonogashira acetylene_acyl halide",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Ullmann condensation",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Goldberg coupling",
        "Suzuki",
        "Stille",
        "decarboxylative_coupling",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "Chan-Lam alcohol",
        "Chan-Lam amine",
        "Chan-Lam etherification",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling

        # Check if this is a reaction node in the late stage (depth <= 1)
        if node["type"] == "reaction" and depth <= 1:
            print(f"Examining reaction at depth {depth}")
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")
                print(f"Number of reactants: {len(reactants)}")

                # For depth 0, verify this is the final product
                if depth == 0:
                    product_mol = Chem.MolFromSmiles(product)
                    final_mol = Chem.MolFromSmiles(final_product_smiles)
                    if product_mol and final_mol:
                        if product_mol.GetNumHeavyAtoms() != final_mol.GetNumHeavyAtoms():
                            print(f"Not the final product step (atom count mismatch)")
                            return

                # Check if this is a known coupling reaction
                is_coupling_reaction = False
                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected coupling reaction: {rxn_type}")
                        is_coupling_reaction = True
                        break

                # Count substantial reactants (excluding small molecules and reagents)
                substantial_reactants = []
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol is None:
                        continue

                    # Consider size, structural complexity, and functional groups
                    heavy_atoms = mol.GetNumHeavyAtoms()
                    ring_count = len(Chem.GetSSSR(mol))

                    # Check for important ring structures
                    has_important_ring = any(
                        checker.check_ring(ring, reactant) for ring in important_rings
                    )

                    # Check for important functional groups
                    has_important_fg = any(checker.check_fg(fg, reactant) for fg in important_fgs)

                    # A substantial fragment has either:
                    # - At least 7 heavy atoms, or
                    # - At least 5 heavy atoms and contains a ring, or
                    # - Contains an important ring structure, or
                    # - Contains an important functional group
                    if (
                        (heavy_atoms >= 7)
                        or (heavy_atoms >= 5 and ring_count > 0)
                        or has_important_ring
                        or has_important_fg
                    ):
                        substantial_reactants.append(reactant)
                        print(
                            f"Substantial reactant: {reactant} (atoms: {heavy_atoms}, rings: {ring_count})"
                        )

                # Check for fragment coupling conditions
                if len(substantial_reactants) >= 2:
                    product_mol = Chem.MolFromSmiles(product)

                    if is_coupling_reaction:
                        print(
                            f"Confirmed late-stage fragment coupling with {len(substantial_reactants)} substantial fragments"
                        )
                        late_coupling = True
                    else:
                        # Check if the product is significantly larger than any single reactant
                        reactant_mols = [
                            Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                        ]
                        if reactant_mols:
                            max_reactant_atoms = max(
                                [mol.GetNumHeavyAtoms() for mol in reactant_mols]
                            )
                            product_atoms = product_mol.GetNumHeavyAtoms() if product_mol else 0

                            if (
                                product_atoms > max_reactant_atoms * 1.3
                            ):  # Product is significantly larger
                                print(
                                    f"Detected fragment joining: product size ({product_atoms}) > 1.3× largest reactant ({max_reactant_atoms})"
                                )

                                # Try to find MCS between product and each substantial reactant
                                # to confirm they're incorporated largely intact
                                intact_fragments = 0
                                for reactant in substantial_reactants:
                                    r_mol = Chem.MolFromSmiles(reactant)
                                    if r_mol and product_mol:
                                        mcs = rdFMCS.FindMCS(
                                            [r_mol, product_mol],
                                            completeRingsOnly=True,
                                            ringMatchesRingOnly=True,
                                        )
                                        # Use both percentage and absolute minimum threshold
                                        min_atoms_threshold = min(r_mol.GetNumHeavyAtoms() * 0.7, 5)
                                        if mcs.numAtoms >= min_atoms_threshold:
                                            intact_fragments += 1
                                            print(
                                                f"Fragment largely intact: {reactant} ({mcs.numAtoms} of {r_mol.GetNumHeavyAtoms()} atoms matched)"
                                            )

                                if intact_fragments >= 2:
                                    print(
                                        f"Found {intact_fragments} substantial fragments largely intact in product"
                                    )
                                    late_coupling = True

                # Check for specific coupling functional groups
                if not late_coupling and len(substantial_reactants) >= 2:
                    # Check for boronic acids/esters (Suzuki coupling)
                    has_boronic = any(
                        checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    )
                    # Check for halides (many coupling reactions)
                    has_halide = any(
                        checker.check_fg("Aromatic halide", r)
                        or checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        or checker.check_fg("Alkenyl halide", r)
                        for r in reactants
                    )
                    # Check for tin compounds (Stille)
                    has_tin = any(checker.check_fg("Tin", r) for r in reactants)
                    # Check for triflates/tosylates (coupling partners)
                    has_leaving_group = any(
                        checker.check_fg("Triflate", r)
                        or checker.check_fg("Tosylate", r)
                        or checker.check_fg("Mesylate", r)
                        for r in reactants
                    )

                    if (has_boronic and (has_halide or has_leaving_group)) or (
                        has_tin and (has_halide or has_leaving_group)
                    ):
                        print(
                            "Detected coupling reaction components (boronic/tin + halide/leaving group)"
                        )
                        late_coupling = True

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    print("Starting traversal")
    dfs_traverse(route)
    print(f"Final result: {late_coupling}")
    return late_coupling
