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
    This function detects if the synthesis follows a linear strategy (as opposed to convergent).
    A linear synthesis typically builds a molecule step by step, with at most one step having
    multiple reactants. A convergent synthesis combines multiple complex fragments.
    """
    multi_reactant_steps = 0
    total_steps = 0

    # Protection/deprotection and simple functional group transformations
    simple_transformation_reactions = [
        "Boc amine protection",
        "Boc amine deprotection",
        "Alcohol protection with silyl ethers",
        "Alcohol deprotection from silyl ethers",
        "Protection of carboxylic acid",
        "Deprotection of carboxylic acid",
        "Ester saponification",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "TMS deprotection from alkyne",
        "Phthalimide deprotection",
        "N-glutarimide deprotection",
        "Cleavage of methoxy ethers to alcohols",
        "Cleavage of alkoxy ethers to alcohols",
        "Ether cleavage to primary alcohol",
        "COOH ethyl deprotection",
        "Tert-butyl deprotection of amine",
        "Oxidation of aldehydes to carboxylic acids",
        "Reduction of ester to primary alcohol",
        "Reduction of ketone to secondary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Nitrile to amide",
        "Reduction of nitrile to amine",
        "Reduction of primary amides to amines",
        "Reduction of secondary amides to amines",
        "Reduction of tertiary amides to amines",
        "Alcohol to chloride",
        "Alcohol to triflate conversion",
        "Carboxylic acid to carboxylate",
        "Ester to carboxylate",
        "Carboxylate to carboxylic acid",
    ]

    # C-C bond forming reactions that are typically used in convergent synthesis
    convergent_reactions = [
        "Suzuki coupling",
        "Stille reaction",
        "Negishi coupling",
        "Heck",
        "Sonogashira",
        "Buchwald-Hartwig",
        "Ullmann-Goldberg",
        "Grignard",
        "Wittig",
        "Diels-Alder",
        "Aldol condensation",
        "Michael addition",
        "Pauson-Khand reaction",
    ]

    def is_simple_reagent(smiles):
        """Check if a reactant is a simple reagent rather than a complex building block"""
        # Simple reagents typically have short SMILES strings
        if len(smiles) < 8:
            return True

        # Check for common simple reagents
        simple_reagents = [
            "O",
            "CO",
            "[OH-]",
            "[H+]",
            "Cl",
            "Br",
            "I",
            "[Na+]",
            "[K+]",
            "CC(=O)O",
            "CCO",
            "CCOCC",
            "CN",
            "C1CCOC1",
            "CC(C)O",
            "CC#N",
            "CS(=O)(=O)O",
            "ClCCl",
            "BrCBr",
            "ICCCI",
            "CC(=O)OC(=O)C",
            "CC(=O)Cl",
            "CC(=O)Br",
            "CS(=O)(=O)Cl",
            "CC(C)(C)OC(=O)Cl",
            "CC(C)(C)OC(=O)O",
            "CF3C(=O)O",
            "C[Si](C)(C)Cl",
            "C[Si](C)(C)OC(C)(C)C",
            "B(O)O",
            "B(OC)OC",
            "B(OCC)OCC",
            "[Pd]",
            "[Pt]",
            "[Rh]",
            "[Ru]",
            "[Cu]",
            "[Zn]",
            "[Mg]",
            "[Li]",
            "[Al]",
        ]
        if smiles in simple_reagents:
            return True

        # Check for common protecting/deprotecting reagents
        protecting_reagents = [
            "CC(C)(C)OC(=O)",
            "CF3C(=O)O",
            "CS(=O)(=O)O",
            "CC(=O)",
            "C[Si](C)(C)",
            "C[Si](C)(C)C",
            "C[Si](CC)(CC)CC",
            "C(=O)OC(=O)",
            "C(=O)Cl",
            "C(=O)Br",
            "S(=O)(=O)Cl",
            "P(=O)(Cl)Cl",
        ]
        for reagent in protecting_reagents:
            if reagent in smiles:
                return True

        # Check for common functional groups that indicate simple reagents
        simple_fgs = [
            "Acyl halide",
            "Triflate",
            "Mesylate",
            "Tosylate",
            "Anhydride",
            "Formaldehyde",
            "Carbon dioxide",
            "Carbon monoxide",
            "Methanol",
        ]
        for fg in simple_fgs:
            if checker.check_fg(fg, smiles):
                return True

        # Check molecular complexity - simple reagents have fewer atoms
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumAtoms() < 10:
            return True

        return False

    def is_simple_transformation(rxn_smiles):
        """Check if the reaction is a simple transformation (protection, deprotection, FG transformation)"""
        # Check predefined simple transformations
        for rxn_type in simple_transformation_reactions:
            if checker.check_reaction(rxn_type, rxn_smiles):
                print(f"Simple transformation detected: {rxn_type}")
                return True

        # Check if it's a functional group transformation
        fg_transformations = [
            "Oxidation",
            "Reduction",
            "Hydrolysis",
            "Esterification",
            "Amidation",
            "Dehydration",
            "Hydration",
            "Halogenation",
            "Dehalogenation",
        ]
        for rxn_pattern in fg_transformations:
            if rxn_pattern.lower() in rxn_smiles.lower():
                print(f"Functional group transformation detected: {rxn_pattern}")
                return True

        # Check if it's a C-C bond forming reaction (typically convergent)
        for rxn_type in convergent_reactions:
            if checker.check_reaction(rxn_type, rxn_smiles):
                print(f"C-C bond forming reaction detected: {rxn_type}")
                return False

        return False

    def dfs_traverse(node, depth=0):
        nonlocal multi_reactant_steps, total_steps

        if node["type"] == "reaction":
            total_steps += 1

            # Extract reactants and reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            rxn_smiles = node["metadata"].get("smiles", "")
            reactants = rsmi.split(">")[0].split(".")

            # Filter out simple reagents
            complex_reactants = [r for r in reactants if not is_simple_reagent(r)]

            # Check if this is a simple transformation step
            is_simple_step = is_simple_transformation(rxn_smiles)

            # Count as multi-reactant step only if:
            # 1. Has multiple complex reactants
            # 2. Not a simple transformation step
            if len(complex_reactants) > 1 and not is_simple_step:
                multi_reactant_steps += 1
                print(
                    f"Multi-reactant step found at depth {depth}: {len(complex_reactants)} complex reactants"
                )
                print(f"Reactants: {complex_reactants}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Calculate percentage of multi-reactant steps
    percent_multi = (multi_reactant_steps / total_steps * 100) if total_steps > 0 else 0

    print(f"Total steps: {total_steps}, Multi-reactant steps: {multi_reactant_steps}")

    # For a truly linear synthesis:
    # - At most one multi-reactant step for short syntheses (≤5 steps)
    # - For longer syntheses, allow a small percentage of multi-reactant steps
    if multi_reactant_steps == 0:
        print("Perfectly linear synthesis detected (no multi-reactant steps)")
        return True
    elif multi_reactant_steps == 1 and total_steps <= 5:
        print("Linear synthesis strategy detected (one multi-reactant step in short synthesis)")
        return True
    elif total_steps > 5 and percent_multi <= 30:  # More lenient threshold
        print(f"Mostly linear synthesis detected ({percent_multi:.1f}% multi-reactant steps)")
        return True
    else:
        print(
            f"Convergent synthesis detected with {multi_reactant_steps} multi-reactant steps ({percent_multi:.1f}%)"
        )
        return False
