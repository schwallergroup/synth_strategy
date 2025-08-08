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
    This function detects if the synthesis involves formation of a heterocyclic system
    in the final step (late-stage heterocycle formation).
    """
    has_late_stage_heterocycle = False

    # List of common heterocycles to check
    heterocycles = [
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
        "indole",
        "quinoline",
        "isoquinoline",
        "thiophene",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    # List of heterocycle-forming reaction types
    heterocycle_forming_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_heterocycle

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            print(f"Examining reaction at depth {depth}: {node['metadata'].get('rsmi', 'No RSMI')}")

            # Check if this is the final reaction (depth 1 in retrosynthetic direction)
            if depth == 1:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                print(f"Final reaction - Reactants: {reactants_part}, Product: {products_part}")

                # First check if this is a known heterocycle-forming reaction
                for rxn_type in heterocycle_forming_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected known heterocycle-forming reaction: {rxn_type}")
                        has_late_stage_heterocycle = True
                        return

                # Parse reactants and product
                reactants = reactants_part.split(".")
                product = products_part

                # Check if any heterocycle is present in the product but not in any reactant
                product_heterocycles = []
                reactant_heterocycles = []

                # Check heterocycles in product
                for cycle in heterocycles:
                    if checker.check_ring(cycle, product):
                        product_heterocycles.append(cycle)
                        print(f"Found {cycle} in product")

                # Check heterocycles in reactants
                for reactant in reactants:
                    reactant_cycles = []
                    for cycle in heterocycles:
                        if checker.check_ring(cycle, reactant):
                            reactant_cycles.append(cycle)
                            reactant_heterocycles.append(cycle)
                    print(f"Found cycles in reactant {reactant}: {reactant_cycles}")

                # Find heterocycles in product that aren't in any reactant
                new_heterocycles = [
                    cycle for cycle in product_heterocycles if cycle not in reactant_heterocycles
                ]

                if new_heterocycles:
                    print(f"Detected late-stage heterocycle formation: {new_heterocycles}")
                    has_late_stage_heterocycle = True
                    return

                # Check for triazole formation specifically (common in click chemistry)
                if checker.check_ring("triazole", product):
                    triazole_in_reactants = any(
                        checker.check_ring("triazole", r) for r in reactants
                    )
                    if not triazole_in_reactants:
                        print("Detected triazole formation")
                        has_late_stage_heterocycle = True
                        return

                # Fallback to ring counting method if no specific heterocycles detected
                try:
                    reactants_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and all(reactants_mols):
                        # Count rings in product
                        product_rings = Chem.GetSSSR(product_mol)

                        # Count total rings in reactants
                        reactant_rings_total = sum(len(Chem.GetSSSR(m)) for m in reactants_mols)

                        print(
                            f"Ring count - Product: {len(product_rings)}, Reactants total: {reactant_rings_total}"
                        )

                        # If product has more rings than reactants combined, a ring was formed
                        if len(product_rings) > reactant_rings_total:
                            # Check if any new ring contains a heteroatom
                            for ring in product_rings:
                                has_heteroatom = False
                                for atom_idx in ring:
                                    atom = product_mol.GetAtomWithIdx(atom_idx)
                                    if atom.GetSymbol() not in ["C", "H"]:
                                        has_heteroatom = True
                                        break

                                if has_heteroatom:
                                    print(
                                        "Detected late-stage heterocycle formation using ring count method"
                                    )
                                    has_late_stage_heterocycle = True
                                    break
                except Exception as e:
                    print(f"Error in ring counting method: {e}")

        # Continue traversing with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_heterocycle
