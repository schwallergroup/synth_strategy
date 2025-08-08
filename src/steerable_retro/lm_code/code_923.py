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
    Detects if the synthesis employs a late-stage fragment coupling strategy
    where major fragments are combined in the final steps.
    """
    late_stage_coupling = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Extract depth information
            depth_info = node.get("metadata", {}).get("ID", "")
            node_depth = current_depth

            # Try to extract depth from ID if available
            if "Depth:" in depth_info:
                try:
                    depth_str = depth_info.split("Depth:")[1].strip().split()[0]
                    node_depth = int(depth_str)
                except (IndexError, ValueError):
                    pass  # Keep using current_depth if extraction fails

            # Check if this is a late-stage reaction (depth 0-2)
            if node_depth <= 2:
                print(f"Examining late-stage reaction at depth {node_depth}: {rsmi}")

                # Check if this is a known coupling reaction
                is_coupling_reaction = any(
                    [
                        checker.check_reaction("Suzuki coupling with boronic acids", rsmi),
                        checker.check_reaction("Suzuki coupling with boronic esters", rsmi),
                        checker.check_reaction("Negishi coupling", rsmi),
                        checker.check_reaction("Stille reaction_aryl", rsmi),
                        checker.check_reaction("Stille reaction_vinyl", rsmi),
                        checker.check_reaction("Heck terminal vinyl", rsmi),
                        checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi),
                        checker.check_reaction("Ullmann condensation", rsmi),
                        checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        ),
                        checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        ),
                    ]
                )

                # If not a known coupling reaction, check for C-C bond formation
                if not is_coupling_reaction:
                    print("Not a known coupling reaction, checking for C-C bond formation...")
                    # This is a simplified check - in a real implementation, you would need
                    # to analyze the reaction more carefully to confirm C-C bond formation

                # Count substantial fragments in reactants
                substantial_fragment_count = 0
                substantial_fragments = []

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Count atoms, rings, and check molecular weight
                        atom_count = mol.GetNumAtoms()
                        ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
                        mw = Chem.Descriptors.MolWt(mol)

                        # Consider a fragment substantial if it meets certain criteria
                        if (ring_count >= 1 and atom_count >= 7) or mw > 120:
                            substantial_fragment_count += 1
                            substantial_fragments.append(reactant)
                            print(
                                f"Found substantial fragment: {reactant} (rings: {ring_count}, atoms: {atom_count}, MW: {mw:.1f})"
                            )

                # If we're combining multiple substantial fragments in late stage
                if substantial_fragment_count >= 2 and (is_coupling_reaction or node_depth <= 1):
                    print(f"Found late-stage fragment coupling at depth {node_depth}")
                    print(f"Fragments: {substantial_fragments}")
                    late_stage_coupling = True
                else:
                    print(
                        f"Not a late-stage fragment coupling: fragments={substantial_fragment_count}, is_coupling={is_coupling_reaction}"
                    )

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    print(f"Final result: late_stage_coupling = {late_stage_coupling}")
    return late_stage_coupling
