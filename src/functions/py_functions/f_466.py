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
    Detects if the synthesis uses a late-stage heterocycle formation strategy,
    specifically forming a heterocyclic ring in the final step.
    """
    # List of heterocyclic rings to check
    heterocycle_rings = [
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

    # List of heterocycle-forming reactions
    heterocycle_forming_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "spiro-chromanone",
        "pyrazole",
        "phthalazinone",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "piperidine_indole",
        "imidazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    ]

    def is_final_reaction(node, depth):
        """Check if this is the final reaction in the synthesis route"""
        # If depth is 0 and node is a reaction, it's the final reaction
        if depth == 0 and node["type"] == "reaction":
            return True

        # If depth is 1, node is a molecule, and it has no parent reaction, it could be the final product
        if depth == 1 and node["type"] == "mol" and len(node.get("children", [])) == 1:
            child = node["children"][0]
            if child["type"] == "reaction":
                return True

        return False

    def check_heterocycle_formation(rsmi, reactants_smiles, product_smiles):
        """Check if the reaction forms a heterocycle"""
        print(f"Checking reaction: {rsmi}")

        # 1. Check if this is a known heterocycle-forming reaction
        for reaction_type in heterocycle_forming_reactions:
            if checker.check_reaction(reaction_type, rsmi):
                print(f"Late-stage heterocycle formation detected: {reaction_type}")
                return True

        # Parse molecules
        try:
            product_mol = Chem.MolFromSmiles(product_smiles)
            if not product_mol:
                print(f"Could not parse product SMILES: {product_smiles}")
                return False

            reactant_mols = []
            for r_smiles in reactants_smiles:
                mol = Chem.MolFromSmiles(r_smiles)
                if mol:
                    reactant_mols.append(mol)
                else:
                    print(f"Could not parse reactant SMILES: {r_smiles}")

            if not reactant_mols:
                print("No valid reactant molecules found")
                return False

            # 2. Check if product contains heterocycles not present in reactants
            for ring_name in heterocycle_rings:
                if checker.check_ring(ring_name, product_smiles):
                    print(f"Product contains {ring_name}")

                    # Check if this ring is new (not in any reactant)
                    ring_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_ring(ring_name, reactant):
                            print(f"Reactant also contains {ring_name}: {reactant}")
                            ring_in_reactants = True
                            break

                    if not ring_in_reactants:
                        print(
                            f"Late-stage heterocycle formation detected: {ring_name} formed in final step"
                        )
                        return True

            # 3. Check if there's an increase in ring count
            product_ring_info = product_mol.GetRingInfo()
            product_ring_count = product_ring_info.NumRings()

            max_reactant_ring_count = max(
                [mol.GetRingInfo().NumRings() for mol in reactant_mols]
            )

            print(
                f"Product ring count: {product_ring_count}, Max reactant ring count: {max_reactant_ring_count}"
            )

            if product_ring_count > max_reactant_ring_count:
                # Check if any of the new rings contains heteroatoms
                for ring_idx in range(product_ring_info.NumRings()):
                    ring_atoms = product_ring_info.AtomRings()[ring_idx]
                    has_heteroatom = False
                    for atom_idx in ring_atoms:
                        atom = product_mol.GetAtomWithIdx(atom_idx)
                        if atom.GetAtomicNum() not in [1, 6]:  # Not H or C
                            has_heteroatom = True
                            break

                    if has_heteroatom:
                        print(
                            "Late-stage heterocycle formation detected: New heterocycle formed in final step"
                        )
                        return True

            return False

        except Exception as e:
            print(f"Error in check_heterocycle_formation: {e}")
            return False

    def dfs_traverse(node, depth=0):
        """Traverse the synthesis route to find late-stage heterocycle formation"""
        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # Check if this is the final or near-final reaction
            if depth <= 1:
                print(f"Examining potential final reaction at depth {depth}")
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if check_heterocycle_formation(rsmi, reactants_smiles, product_smiles):
                    return True

        # Continue traversing
        for child in node.get("children", []):
            if dfs_traverse(child, depth + 1):
                return True

        return False

    # Start traversal
    result = dfs_traverse(route)
    print(f"Final result: {result}")
    return result
