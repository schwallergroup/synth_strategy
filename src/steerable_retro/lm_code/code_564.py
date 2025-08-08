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
    This function detects a linear synthesis strategy with heterocycle formation
    and late-stage fragment coupling.
    """
    # Track key events
    heterocycle_formation = False
    late_coupling = False
    linear_structure = True

    # Track branching factor
    max_children = 0

    # List of common heterocycles to check
    heterocycles = [
        "furan",
        "pyrrole",
        "pyridine",
        "thiophene",
        "oxazole",
        "thiazole",
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "morpholine",
        "piperidine",
        "piperazine",
        "isoxazole",
        "pyrimidine",
        "pyrazine",
        "thiazole",
        "oxadiazole",
        "thiadiazole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
    ]

    # List of common coupling reactions
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Heck terminal vinyl",
        "Ullmann condensation",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Stille reaction_vinyl",
        "Suzuki coupling with sulfonic esters",
        "Heck_non-terminal_vinyl",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Goldberg coupling",
        "Catellani reaction ortho",
        "Catellani reaction para",
        "Aryllithium cross-coupling",
        "{Suzuki}",
        "{Heck_terminal_vinyl}",
        "{Heck_non-terminal_vinyl}",
        "{Stille}",
        "{N-arylation_heterocycles}",
        "{Buchwald-Hartwig}",
        "{decarboxylative_coupling}",
    ]

    # Additional fragment coupling reactions
    fragment_coupling_reactions = [
        "Ugi reaction",
        "Petasis reaction",
        "A3 coupling",
        "Aldol condensation",
        "Wittig reaction",
        "Julia Olefination",
        "Michael addition",
        "Diels-Alder",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "{reductive amination}",
        "{Wittig}",
        "{Mitsunobu_imide}",
        "{Mitsunobu_phenole}",
        "{Mitsunobu_sulfonamide}",
        "{Mitsunobu_tetrazole_1}",
        "{Mitsunobu_tetrazole_2}",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation, late_coupling, linear_structure, max_children

        if node["type"] == "reaction":
            try:
                # Get reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Count reactants to check for convergent steps
                if len(reactants_smiles) > max_children:
                    max_children = len(reactants_smiles)

                # Check for heterocycle formation (middle of synthesis)
                if 0 < depth < 4 and not heterocycle_formation:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    reactant_mols = [mol for mol in reactant_mols if mol is not None]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol:
                        # Count rings in reactants and product
                        reactant_rings = sum(
                            [mol.GetRingInfo().NumRings() for mol in reactant_mols if mol]
                        )
                        product_rings = product_mol.GetRingInfo().NumRings()

                        # Check if a new heterocycle was formed
                        if product_rings > reactant_rings:
                            # Verify it's a heterocycle by checking against known heterocycles
                            for heterocycle in heterocycles:
                                if checker.check_ring(heterocycle, product_smiles):
                                    # Check if this heterocycle wasn't in any of the reactants
                                    if not any(
                                        checker.check_ring(heterocycle, r) for r in reactants_smiles
                                    ):
                                        print(f"Detected {heterocycle} formation at depth {depth}")
                                        heterocycle_formation = True
                                        break

                # Check for late-stage coupling (final step or one step before)
                if depth <= 1:
                    if len(reactants_smiles) >= 2:
                        # Verify it's a coupling reaction
                        for coupling in coupling_reactions:
                            if checker.check_reaction(coupling, rsmi):
                                print(f"Detected late-stage {coupling} at depth {depth}")
                                late_coupling = True
                                break

                        # If no specific coupling reaction was found, check for other fragment coupling reactions
                        if not late_coupling:
                            for rxn_name in fragment_coupling_reactions:
                                if checker.check_reaction(rxn_name, rsmi):
                                    print(f"Detected late-stage {rxn_name} at depth {depth}")
                                    late_coupling = True
                                    break

                        # Final fallback: any reaction with multiple reactants at depth 0-1
                        if not late_coupling:
                            print(
                                f"Detected late-stage fragment coupling (multiple reactants) at depth {depth}"
                            )
                            late_coupling = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis has limited branching
    if max_children > 2:
        linear_structure = False

    print(
        f"Strategy analysis: heterocycle_formation={heterocycle_formation}, late_coupling={late_coupling}, linear_structure={linear_structure}"
    )

    # Return True if all key elements of the strategy are present
    return heterocycle_formation and late_coupling and linear_structure
