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
    This function detects if the synthesis involves late-stage functional group
    modifications (in the final steps) while maintaining the core scaffold.
    """
    fg_modifications_at_depth = {}
    max_depth = -1

    # List of functional groups to check
    fg_list = [
        "Amine",
        "Nitro group",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Alcohol",
        "Aldehyde",
        "Ketone",
        "Nitrile",
        "Azide",
        "Alkyne",
        "Alkene",
        "Ether",
        "Thiol",
        "Sulfide",
        "Sulfone",
        "Sulfoxide",
        "Phosphate ester",
        "Boronic acid",
        "Boronic ester",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            product_mol = Chem.MolFromSmiles(product_part)
            if not product_mol:
                print(f"Failed to parse product SMILES: {product_part}")
                return

            # Split reactants
            reactant_smiles = reactants_part.split(".")
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactant_smiles if r]
            reactant_mols = [m for m in reactant_mols if m]  # Filter out None values

            if not reactant_mols:
                print(f"Failed to parse any reactant SMILES from: {reactants_part}")
                return

            # Check for scaffold preservation using MCS
            for reactant_mol in reactant_mols:
                mcs = rdFMCS.FindMCS(
                    [reactant_mol, product_mol],
                    completeRingsOnly=True,
                    ringMatchesRingOnly=True,
                    matchValences=True,
                )

                if mcs.numAtoms > 0:
                    # Calculate scaffold similarity
                    scaffold_ratio = mcs.numAtoms / product_mol.GetNumAtoms()

                    # If significant scaffold is preserved (>60%), check for FG modifications
                    if scaffold_ratio > 0.6:
                        # Check for functional group changes
                        fg_modified = False
                        for fg in fg_list:
                            reactant_has_fg = any(
                                checker.check_fg(fg, Chem.MolToSmiles(r)) for r in reactant_mols
                            )
                            product_has_fg = checker.check_fg(fg, product_part)

                            if reactant_has_fg != product_has_fg:
                                if depth not in fg_modifications_at_depth:
                                    fg_modifications_at_depth[depth] = []
                                fg_modifications_at_depth[depth].append(fg)
                                print(f"Found {fg} modification at depth {depth}: {rsmi}")
                                fg_modified = True

                        # Also check for specific reaction types that indicate FG modifications
                        rxn_types = [
                            "Oxidation of aldehydes to carboxylic acids",
                            "Reduction of aldehydes and ketones to alcohols",
                            "Reduction of nitrile to amine",
                            "Esterification of Carboxylic Acids",
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        ]

                        for rxn_type in rxn_types:
                            if checker.check_reaction(rxn_type, rsmi):
                                if depth not in fg_modifications_at_depth:
                                    fg_modifications_at_depth[depth] = []
                                fg_modifications_at_depth[depth].append(f"Reaction: {rxn_type}")
                                print(f"Found reaction {rxn_type} at depth {depth}: {rsmi}")
                                fg_modified = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if functional group modifications occur in the late stage (first half of synthesis)
    late_stage_modifications = False
    if max_depth > 0:
        late_stage_threshold = max_depth // 2
        # In retrosynthetic analysis, lower depths are later stages
        late_stage_modifications = any(
            depth <= late_stage_threshold for depth in fg_modifications_at_depth.keys()
        )
        print(f"Max depth: {max_depth}, Late stage threshold: {late_stage_threshold}")
        print(f"FG modifications at depths: {list(fg_modifications_at_depth.keys())}")
        print(f"Late stage modifications: {late_stage_modifications}")

    return late_stage_modifications
