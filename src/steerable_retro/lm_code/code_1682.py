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
    This function detects if the synthetic route involves sequential elaboration
    of heterocyclic systems through multiple functionalization steps.
    """
    heterocycle_steps = 0
    sequential_steps = 0
    current_sequence = 0

    # List of heterocyclic rings to check
    heterocyclic_rings = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "trioxane",
        "dioxepane",
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
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "dithiane",
        "dithiolane",
        "benzothiophene",
        "oxathiolane",
        "dioxathiolane",
        "thiazolidine",
        "oxazolidine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "pteridin",
        "phenothiazine",
        "phenoxazine",
        "dibenzofuran",
        "dibenzothiophene",
        "xanthene",
        "thioxanthene",
        "pyrroline",
        "pyrrolidone",
        "imidazolidine",
        "porphyrin",
        "indazole",
        "benzotriazole",
    ]

    # Reactions that typically modify heterocycles
    heterocycle_modifying_reactions = [
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "Aromatic halogenation",
        "Aromatic nitration",
        "N-alkylation",
        "Aromatic substitution",
        "Minisci",
        "Directed ortho metalation",
        "Suzuki coupling",
        "Buchwald-Hartwig",
        "Heck",
        "Sonogashira",
        "Negishi coupling",
        "Stille reaction",
        "Chan-Lam",
    ]

    def has_heterocycle(mol_smiles):
        """Check if molecule contains any heterocyclic ring"""
        for ring in heterocyclic_rings:
            if checker.check_ring(ring, mol_smiles):
                return True
        return False

    def is_heterocycle_modifying_reaction(rxn_smiles):
        """Check if reaction is likely modifying a heterocycle"""
        for rxn_type in heterocycle_modifying_reactions:
            if checker.check_reaction(rxn_type, rxn_smiles):
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_steps, sequential_steps, current_sequence

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                try:
                    # Extract reactants and product
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")

                    # Check if product contains heterocycle
                    product_has_heterocycle = has_heterocycle(product_part)

                    # Check if any reactant contains heterocycle
                    reactant_has_heterocycle = any(has_heterocycle(r) for r in reactants)

                    # Check if this is a heterocycle modifying reaction
                    is_modifying_reaction = is_heterocycle_modifying_reaction(rsmi)

                    # Determine if this is a heterocycle elaboration step
                    if product_has_heterocycle and (
                        reactant_has_heterocycle or is_modifying_reaction
                    ):
                        heterocycle_steps += 1
                        current_sequence += 1
                        print(f"Heterocycle modification at depth {depth}: {rsmi}")
                    else:
                        # Reset sequence counter if we encounter a non-heterocycle step
                        sequential_steps = max(sequential_steps, current_sequence)
                        current_sequence = 0
                        print(f"Non-heterocycle step at depth {depth}: {rsmi}")
                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")
                    # Don't reset sequence on error

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Update sequential_steps one last time in case the sequence continues to the end
    sequential_steps = max(sequential_steps, current_sequence)

    has_sequential_elaboration = sequential_steps >= 3
    print(f"Heterocycle steps: {heterocycle_steps}")
    print(f"Longest sequential elaboration: {sequential_steps}")
    print(f"Sequential heterocycle elaboration: {has_sequential_elaboration}")

    return has_sequential_elaboration
