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
    Detects synthesis of thiophene-fused pyrimidine systems where a nucleophilic
    aromatic substitution (SNAr) is used to introduce a nitrogen heterocycle.
    """
    # Track if we found SNAr on pyrimidine
    found_snar = False
    # Track if thiophene-pyrimidine system is present in final product
    has_thiophene_pyrimidine = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar, has_thiophene_pyrimidine

        # Check if final product has thiophene-pyrimidine system
        if node["type"] == "mol" and depth == 0:
            mol_smiles = node["smiles"]
            print(f"Checking final product: {mol_smiles}")
            # Check if both thiophene and pyrimidine are present in the molecule
            if checker.check_ring("thiophene", mol_smiles) and checker.check_ring(
                "pyrimidine", mol_smiles
            ):
                print("Both thiophene and pyrimidine rings found in final product")
                # Get atom indices for both rings
                try:
                    from rdkit import Chem

                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        # Get ring atom indices
                        thiophene_indices = checker.get_ring_atom_indices("thiophene", mol_smiles)
                        pyrimidine_indices = checker.get_ring_atom_indices("pyrimidine", mol_smiles)

                        # Extract atom indices from potentially nested structures
                        thiophene_atoms = set()
                        pyrimidine_atoms = set()

                        # Process thiophene indices
                        if isinstance(thiophene_indices, list):
                            for item in thiophene_indices:
                                if isinstance(item, tuple):
                                    for subitem in item:
                                        if isinstance(subitem, int):
                                            thiophene_atoms.add(subitem)
                                elif isinstance(item, int):
                                    thiophene_atoms.add(item)
                        elif isinstance(thiophene_indices, int):
                            thiophene_atoms.add(thiophene_indices)

                        # Process pyrimidine indices
                        if isinstance(pyrimidine_indices, list):
                            for item in pyrimidine_indices:
                                if isinstance(item, tuple):
                                    for subitem in item:
                                        if isinstance(subitem, int):
                                            pyrimidine_atoms.add(subitem)
                                elif isinstance(item, int):
                                    pyrimidine_atoms.add(item)
                        elif isinstance(pyrimidine_indices, int):
                            pyrimidine_atoms.add(pyrimidine_indices)

                        print(f"Thiophene atoms: {thiophene_atoms}")
                        print(f"Pyrimidine atoms: {pyrimidine_atoms}")

                        # Check if any thiophene atom is bonded to any pyrimidine atom
                        for t_atom in thiophene_atoms:
                            for p_atom in pyrimidine_atoms:
                                bond = mol.GetBondBetweenAtoms(int(t_atom), int(p_atom))
                                if bond is not None:
                                    has_thiophene_pyrimidine = True
                                    print(
                                        f"Found thiophene-pyrimidine system in final product: bond between atoms {t_atom} and {p_atom}"
                                    )
                                    break
                            if has_thiophene_pyrimidine:
                                break
                except Exception as e:
                    print(f"Error processing ring indices: {e}")

        # Check reactions for SNAr
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is a nucleophilic aromatic substitution reaction
                snar_reactions = [
                    "heteroaromatic_nuc_sub",
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                ]
                is_snar = any(checker.check_reaction(rxn, rsmi) for rxn in snar_reactions)

                # If not directly identified as SNAr, check for characteristic patterns
                if not is_snar:
                    # Extract reactants and product
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product has pyrimidine
                    if checker.check_ring("pyrimidine", product):
                        # Check for halide on pyrimidine in reactants and amine nucleophile
                        has_halopyrimidine = False
                        has_amine_nucleophile = False

                        for reactant in reactants:
                            if checker.check_ring("pyrimidine", reactant) and checker.check_fg(
                                "Aromatic halide", reactant
                            ):
                                has_halopyrimidine = True
                                print(f"Found halopyrimidine reactant: {reactant}")

                            # Check for amine nucleophiles
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Tertiary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                                or checker.check_ring("morpholine", reactant)
                                or checker.check_ring("piperidine", reactant)
                                or checker.check_ring("piperazine", reactant)
                            ):
                                has_amine_nucleophile = True
                                print(f"Found amine nucleophile reactant: {reactant}")

                        # If we have both halopyrimidine and amine nucleophile, this is likely SNAr
                        if has_halopyrimidine and has_amine_nucleophile:
                            # Check if the product has an amine connected to pyrimidine where halide was
                            from rdkit import Chem

                            try:
                                # Check if the reaction involves thiophene
                                has_thiophene = any(
                                    checker.check_ring("thiophene", r) for r in reactants
                                ) and checker.check_ring("thiophene", product)

                                if has_thiophene:
                                    print(
                                        f"Identified SNAr reaction with thiophene involvement at depth {depth}"
                                    )
                                    found_snar = True
                            except Exception as e:
                                print(f"Error analyzing reaction: {e}")

                # If already identified as SNAr, check if it involves thiophene and pyrimidine
                elif is_snar:
                    print(f"Found SNAr reaction: {rsmi}")
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product has pyrimidine and thiophene
                    if checker.check_ring("pyrimidine", product) and checker.check_ring(
                        "thiophene", product
                    ):
                        # Check if any reactant has pyrimidine with halide
                        has_halopyrimidine = any(
                            checker.check_ring("pyrimidine", r)
                            and checker.check_fg("Aromatic halide", r)
                            for r in reactants
                        )

                        # Check if any reactant has thiophene
                        has_thiophene = any(checker.check_ring("thiophene", r) for r in reactants)

                        if has_halopyrimidine and has_thiophene:
                            print(
                                f"SNAr reaction involves thiophene and pyrimidine at depth {depth}"
                            )
                            found_snar = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Final result: found_snar={found_snar}, has_thiophene_pyrimidine={has_thiophene_pyrimidine}"
    )
    return found_snar and has_thiophene_pyrimidine
