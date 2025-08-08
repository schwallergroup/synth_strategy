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
    Detects if the synthetic route involves a convergent synthesis strategy
    where two complex fragments are joined in a late-stage coupling reaction.
    """
    convergent_synthesis = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Focus on late-stage reactions (low depth in retrosynthetic tree)
            if depth <= 3:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Only consider reactions with multiple reactants
                if len(reactants) >= 2:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                        # Check if this is a coupling reaction
                        is_coupling = False
                        coupling_reactions = [
                            # C-C coupling reactions
                            "Suzuki coupling with boronic acids",
                            "Suzuki coupling with boronic esters",
                            "Suzuki coupling with sulfonic esters",
                            "Negishi coupling",
                            "Stille reaction_vinyl",
                            "Stille reaction_aryl",
                            "Stille reaction_benzyl",
                            "Stille reaction_allyl",
                            "Heck terminal vinyl",
                            "Sonogashira acetylene_aryl halide",
                            "Sonogashira alkyne_aryl halide",
                            "Kumada cross-coupling",
                            "Hiyama-Denmark Coupling",
                            "Aryllithium cross-coupling",
                            "decarboxylative_coupling",
                            # C-N coupling reactions
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                            "N-arylation_heterocycles",
                            "Goldberg coupling aryl amine-aryl chloride",
                            "Goldberg coupling",
                            # General coupling categories
                            "Suzuki",
                            "Negishi",
                            "Stille",
                            "Heck",
                            "Sonogashira",
                            "Buchwald-Hartwig",
                            "Ullmann",
                            "Kumada",
                        ]

                        for rxn_type in coupling_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                is_coupling = True
                                print(f"Found coupling reaction: {rxn_type}")
                                break

                        # Additional check for general coupling if no specific reaction detected
                        if not is_coupling:
                            # Check for common coupling functional groups
                            has_coupling_fg = False
                            for r in reactants:
                                if r and (
                                    checker.check_fg("Aromatic halide", r)
                                    or checker.check_fg("Boronic acid", r)
                                    or checker.check_fg("Boronic ester", r)
                                ):
                                    has_coupling_fg = True
                                    print(f"Found coupling functional group in: {r}")
                                    break

                            if has_coupling_fg:
                                is_coupling = True
                                print("Detected coupling based on functional groups")

                        # Check complexity of reactants
                        complex_fragments = []
                        for i, r_mol in enumerate(reactant_mols):
                            if r_mol is None:
                                continue

                            # Consider fragments with >8 atoms as potentially complex
                            if r_mol.GetNumAtoms() > 8:
                                # Check for ring structures to confirm complexity
                                has_ring = False
                                ring_info = r_mol.GetRingInfo()
                                if ring_info.NumRings() > 0:
                                    has_ring = True

                                # Check for functional groups to confirm complexity
                                has_fg = False
                                for fg in [
                                    "Aromatic halide",
                                    "Boronic acid",
                                    "Boronic ester",
                                    "Carboxylic acid",
                                    "Ester",
                                    "Amide",
                                    "Amine",
                                    "Aromatic alcohol",
                                    "Nitrile",
                                    "Alkyne",
                                    "Alkene",
                                ]:
                                    if checker.check_fg(fg, Chem.MolToSmiles(r_mol)):
                                        has_fg = True
                                        break

                                if has_ring or has_fg:
                                    complex_fragments.append((i, r_mol))
                                    print(f"Found complex fragment: {Chem.MolToSmiles(r_mol)}")

                        # If we have at least 2 complex fragments being joined in a coupling reaction
                        if len(complex_fragments) >= 2 and is_coupling:
                            # We've found what we're looking for - a convergent synthesis
                            convergent_synthesis = True
                            print(f"Found convergent synthesis at depth {depth}: {rsmi}")
                    except Exception as e:
                        print(f"Error processing reaction SMILES: {rsmi}, Error: {e}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return convergent_synthesis
