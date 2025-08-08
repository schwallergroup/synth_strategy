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
    This function detects a heterocycle formation strategy involving nitrogen and sulfur.
    """
    has_heterocycle_formation = False

    # List of N-S containing heterocycles to check
    ns_heterocycles = [
        "thiazole",
        "benzothiazole",
        "isothiazole",
        "thiadiazole",
        "phenothiazine",
        "thiazolidine",
    ]

    # List of relevant heterocycle formation reactions
    heterocycle_reactions = [
        "benzothiazole formation from aldehyde",
        "benzothiazole formation from acyl halide",
        "benzothiazole formation from ester/carboxylic acid",
        "thiazole",
        "{benzothiazole}",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction SMILES at depth {depth}: {rsmi}")

            try:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a known heterocycle formation reaction
                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found heterocycle formation reaction: {reaction_type}")
                        has_heterocycle_formation = True
                        return

                # Check if product contains an N-S heterocycle that's not in reactants
                product_has_ns_heterocycle = False
                reactants_have_same_heterocycle = False

                # Check product for N-S heterocycles
                for ring_name in ns_heterocycles:
                    if checker.check_ring(ring_name, product_smiles):
                        print(f"Product contains {ring_name}")
                        product_has_ns_heterocycle = True

                        # Check if any reactant has the same heterocycle
                        for reactant_smiles in reactants_smiles:
                            if checker.check_ring(ring_name, reactant_smiles):
                                print(f"Reactant also contains {ring_name}")
                                reactants_have_same_heterocycle = True
                                break

                        if not reactants_have_same_heterocycle:
                            print(f"Found heterocycle formation: {ring_name} formed in product")
                            has_heterocycle_formation = True
                            return

                # If no specific ring was found, try a more general approach
                if not product_has_ns_heterocycle:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        # Check if product has any ring with both N and S
                        ring_info = product_mol.GetRingInfo()
                        for ring_atoms in ring_info.AtomRings():
                            has_nitrogen = False
                            has_sulfur = False

                            for atom_idx in ring_atoms:
                                atom = product_mol.GetAtomWithIdx(atom_idx)
                                if atom.GetSymbol() == "N":
                                    has_nitrogen = True
                                elif atom.GetSymbol() == "S":
                                    has_sulfur = True

                            if has_nitrogen and has_sulfur:
                                print("Found ring with both N and S in product")

                                # Check if reactants don't have this ring pattern
                                ring_in_reactants = False
                                for reactant_smiles in reactants_smiles:
                                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                                    if reactant_mol:
                                        reactant_ring_info = reactant_mol.GetRingInfo()
                                        for reactant_ring in reactant_ring_info.AtomRings():
                                            r_has_nitrogen = False
                                            r_has_sulfur = False

                                            for atom_idx in reactant_ring:
                                                atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                                if atom.GetSymbol() == "N":
                                                    r_has_nitrogen = True
                                                elif atom.GetSymbol() == "S":
                                                    r_has_sulfur = True

                                            if r_has_nitrogen and r_has_sulfur:
                                                ring_in_reactants = True
                                                break

                                if not ring_in_reactants:
                                    print("Found heterocycle formation with N and S")
                                    has_heterocycle_formation = True
                                    return
            except Exception as e:
                print(f"Error processing reaction SMILES: {str(e)}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_heterocycle_formation
