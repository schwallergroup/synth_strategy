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
    Detects if the synthesis follows a convergent approach where two complex fragments
    are prepared separately and combined in a late-stage coupling.
    """
    # Track fragment complexity at coupling step
    found_convergent_synthesis = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_synthesis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a late-stage reaction (depth <= 3)
            if depth <= 3:
                print(f"Checking late-stage reaction at depth {depth}: {rsmi}")

                # Check if this is a coupling reaction
                is_coupling = False
                coupling_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Negishi coupling",
                    "Stille reaction_aryl",
                    "Sonogashira alkyne_aryl halide",
                    "Heck terminal vinyl",
                    "Buchwald-Hartwig",
                    "Ullmann condensation",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                ]

                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found coupling reaction: {rxn_type}")
                        is_coupling = True
                        break

                # If it's a coupling reaction, check for complex fragments
                if is_coupling or len(reactants) >= 2:
                    complex_reactants = 0
                    for r in reactants:
                        try:
                            mol = Chem.MolFromSmiles(r)
                            if mol:
                                # Define complexity criteria
                                atom_count = mol.GetNumAtoms()
                                heavy_atom_count = mol.GetNumHeavyAtoms()
                                ring_count = len(Chem.GetSSSR(mol))

                                # Check for heterocycles
                                has_heterocycle = False
                                heterocycles = [
                                    "pyrazole",
                                    "pyridine",
                                    "pyrimidine",
                                    "imidazole",
                                    "oxazole",
                                    "thiazole",
                                    "furan",
                                    "thiophene",
                                    "pyrrole",
                                ]

                                for ring in heterocycles:
                                    if checker.check_ring(ring, r):
                                        print(f"Found heterocycle {ring} in reactant: {r}")
                                        has_heterocycle = True
                                        break

                                # Check for functional groups that indicate complexity
                                has_complex_fg = False
                                complex_fgs = [
                                    "Ester",
                                    "Amide",
                                    "Nitrile",
                                    "Nitro group",
                                    "Sulfonamide",
                                ]

                                for fg in complex_fgs:
                                    if checker.check_fg(fg, r):
                                        print(
                                            f"Found complex functional group {fg} in reactant: {r}"
                                        )
                                        has_complex_fg = True
                                        break

                                # Define what makes a fragment "complex"
                                is_complex = (
                                    (heavy_atom_count >= 10 and ring_count >= 1)
                                    or (ring_count >= 2)
                                    or (has_heterocycle and heavy_atom_count >= 8)
                                    or (has_complex_fg and heavy_atom_count >= 8)
                                )

                                if is_complex:
                                    complex_reactants += 1
                                    print(
                                        f"Found complex reactant: {r} (atoms: {heavy_atom_count}, rings: {ring_count})"
                                    )
                        except Exception as e:
                            print(f"Error processing reactant {r}: {e}")
                            continue

                    if complex_reactants >= 2:
                        print(
                            f"Found convergent synthesis with {complex_reactants} complex fragments"
                        )
                        found_convergent_synthesis = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_convergent_synthesis
