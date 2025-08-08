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
    Detects synthesis strategy involving heterocycle formation.
    """
    heterocycle_formation = False

    # List of common heterocyclic rings to check
    heterocyclic_rings = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "furan",
        "thiophene",
        "pyran",
        "thiopyran",
        "oxirane",
        "aziridine",
        "oxetane",
        "azetidine",
        "thietane",
        "pyrrolidine",
        "tetrahydrofuran",
        "thiolane",
        "piperidine",
        "tetrahydropyran",
        "thiane",
        "morpholine",
        "thiomorpholine",
        "piperazine",
        "oxazoline",
        "thiazolidine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "triazole",
        "tetrazole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
    ]

    # List of common heterocycle formation reactions
    heterocycle_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "Fischer indole",
        "indole",
        "oxadiazole",
        "imidazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reaction at depth {depth}: {rsmi}")

                # Check if this is a known heterocycle formation reaction
                for reaction_name in heterocycle_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(
                            f"Found heterocycle formation reaction: {reaction_name} at depth {depth}"
                        )
                        heterocycle_formation = True
                        return

                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if all(r is not None for r in reactants) and product is not None:
                    # Count rings in reactants and product
                    reactant_ring_count = sum(r.GetRingInfo().NumRings() for r in reactants)
                    product_ring_count = product.GetRingInfo().NumRings()

                    print(
                        f"Reactant ring count: {reactant_ring_count}, Product ring count: {product_ring_count}"
                    )

                    # Check if a ring was formed or transformed
                    if product_ring_count >= reactant_ring_count:
                        # Check for heterocyclic rings in the product that weren't in the reactants
                        for ring_name in heterocyclic_rings:
                            if checker.check_ring(ring_name, product_smiles):
                                # Check if this ring was already present in any reactant
                                if not any(
                                    checker.check_ring(ring_name, r) for r in reactants_smiles
                                ):
                                    print(f"Found new {ring_name} formation at depth {depth}")
                                    heterocycle_formation = True
                                    return

                        # Check for general heterocyclic structures if specific rings weren't found
                        # Check for nitrogen-containing heterocycles
                        if product.HasSubstructMatch(Chem.MolFromSmarts("[N,n]@[*]")):
                            # Check if this N-containing ring is new
                            n_in_ring_product = product.HasSubstructMatch(
                                Chem.MolFromSmarts("[N,n]@[*]")
                            )
                            n_in_ring_reactants = any(
                                r.HasSubstructMatch(Chem.MolFromSmarts("[N,n]@[*]"))
                                for r in reactants
                                if r is not None
                            )

                            if n_in_ring_product and not n_in_ring_reactants:
                                print(f"Found nitrogen heterocycle formation at depth {depth}")
                                heterocycle_formation = True
                                return

                        # Check for oxygen-containing heterocycles
                        if product.HasSubstructMatch(Chem.MolFromSmarts("[O,o]@[*]")):
                            # Check if this O-containing ring is new
                            o_in_ring_product = product.HasSubstructMatch(
                                Chem.MolFromSmarts("[O,o]@[*]")
                            )
                            o_in_ring_reactants = any(
                                r.HasSubstructMatch(Chem.MolFromSmarts("[O,o]@[*]"))
                                for r in reactants
                                if r is not None
                            )

                            if o_in_ring_product and not o_in_ring_reactants:
                                print(f"Found oxygen heterocycle formation at depth {depth}")
                                heterocycle_formation = True
                                return

                        # Check for sulfur-containing heterocycles
                        if product.HasSubstructMatch(Chem.MolFromSmarts("[S,s]@[*]")):
                            # Check if this S-containing ring is new
                            s_in_ring_product = product.HasSubstructMatch(
                                Chem.MolFromSmarts("[S,s]@[*]")
                            )
                            s_in_ring_reactants = any(
                                r.HasSubstructMatch(Chem.MolFromSmarts("[S,s]@[*]"))
                                for r in reactants
                                if r is not None
                            )

                            if s_in_ring_product and not s_in_ring_reactants:
                                print(f"Found sulfur heterocycle formation at depth {depth}")
                                heterocycle_formation = True
                                return

                # Special case: Check for cyclization reactions that form heterocycles
                if product_ring_count > reactant_ring_count:
                    print(f"Ring formation detected at depth {depth}")

                    # Check if the product contains a heterocycle
                    has_n = product.HasSubstructMatch(Chem.MolFromSmarts("[N,n]@[*]"))
                    has_o = product.HasSubstructMatch(Chem.MolFromSmarts("[O,o]@[*]"))
                    has_s = product.HasSubstructMatch(Chem.MolFromSmarts("[S,s]@[*]"))

                    if has_n or has_o or has_s:
                        print(f"Heterocycle formation detected at depth {depth}")
                        heterocycle_formation = True
                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return heterocycle_formation
