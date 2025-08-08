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
    This function detects late-stage introduction of aromatic substituents.
    """
    late_stage_aromatic = False
    max_depth = 0

    # List of common aromatic rings to check
    aromatic_rings = [
        "benzene",
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "pyrazole",
        "oxazole",
        "thiazole",
        "indole",
        "naphthalene",
    ]

    # List of common aromatic substitution reactions
    aromatic_subst_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Heck terminal vinyl",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Sonogashira alkyne_aryl halide",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic fluorination",
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "Aromatic nitration with HNO3",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Esterification of Carboxylic Acids",
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    ]

    # List of functional groups commonly added to aromatic rings
    aromatic_fgs = [
        "Aromatic halide",
        "Nitro group",
        "Aniline",
        "Phenol",
        "Triflate",
        "Tosylate",
        "Mesylate",
        "Nitrile",
        "Carboxylic acid",
        "Ester",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_aromatic, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is in the late stage of the synthesis (depth <= 2)
            if depth <= 2:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for known aromatic substitution reactions
                for reaction_type in aromatic_subst_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Late-stage aromatic substitution reaction detected: {reaction_type}"
                        )

                        # Verify that the reaction involves an aromatic ring
                        for ring in aromatic_rings:
                            if checker.check_ring(ring, product):
                                print(f"Confirmed: reaction modifies a {ring} ring")
                                late_stage_aromatic = True
                                break

                # Check for changes in aromatic functional groups
                for fg in aromatic_fgs:
                    if checker.check_fg(fg, product):
                        # Check if this functional group was introduced or modified
                        fg_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg(fg, reactant):
                                fg_in_reactants = True
                                break

                        if not fg_in_reactants:
                            # Check if the product has an aromatic ring
                            for ring in aromatic_rings:
                                if checker.check_ring(ring, product):
                                    print(
                                        f"Late-stage aromatic substitution detected: {fg} introduced to {ring} at depth {depth}"
                                    )
                                    late_stage_aromatic = True
                                    break

                # Check for modification of existing aromatic rings
                for ring in aromatic_rings:
                    if checker.check_ring(ring, product):
                        # Check if the ring exists in any reactant
                        for reactant in reactants:
                            if checker.check_ring(ring, reactant):
                                # The ring exists in both reactant and product
                                # Now check if any functional group was added/modified on this ring
                                for fg in aromatic_fgs:
                                    if checker.check_fg(fg, product) and not checker.check_fg(
                                        fg, reactant
                                    ):
                                        print(
                                            f"Late-stage aromatic substitution detected: {fg} added to existing {ring} at depth {depth}"
                                        )
                                        late_stage_aromatic = True
                                    elif not checker.check_fg(fg, product) and checker.check_fg(
                                        fg, reactant
                                    ):
                                        print(
                                            f"Late-stage aromatic substitution detected: {fg} removed from {ring} at depth {depth}"
                                        )
                                        late_stage_aromatic = True

                # Special case: Check for ester hydrolysis to carboxylic acid on aromatic rings
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    for ring in aromatic_rings:
                        if checker.check_ring(ring, product) and checker.check_fg(
                            "Carboxylic acid", product
                        ):
                            for reactant in reactants:
                                if checker.check_ring(ring, reactant) and checker.check_fg(
                                    "Ester", reactant
                                ):
                                    print(
                                        f"Late-stage aromatic substitution detected: Ester hydrolysis on {ring} at depth {depth}"
                                    )
                                    late_stage_aromatic = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Maximum depth found: {max_depth}")
    print(f"Late-stage aromatic substitution detected: {late_stage_aromatic}")
    return late_stage_aromatic
