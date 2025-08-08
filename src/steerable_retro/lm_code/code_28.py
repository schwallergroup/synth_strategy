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
    This function detects a synthetic strategy involving heterocyclic ring opening.
    It looks for reactions where a heterocyclic ring in the reactant is opened in the product.
    """
    ring_opening_detected = False

    # List of common heterocyclic rings to check
    heterocyclic_rings = [
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
        "aziridine",
        "azetidine",
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # List of common ring-opening reaction types
    ring_opening_reactions = [
        "Ring opening of epoxide with amine",
        "Acetal hydrolysis to diol",
        "Acetal hydrolysis to aldehyde",
        "Ketal hydrolysis to ketone",
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        "Oxirane functionalization with ketones",
        "Oxirane functionalization with alkyl iodide",
        "Oxirane functionalization with silyl chloride",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal ring_opening_detected

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # First, check if this is a known ring-opening reaction type
                for rxn_type in ring_opening_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected known ring-opening reaction: {rxn_type}")
                        ring_opening_detected = True
                        return

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if not all(reactant_mols) or not product_mol:
                    print("Warning: Could not parse some molecules")
                    return

                # Count rings in reactants and product
                reactant_ring_counts = [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                product_ring_count = product_mol.GetRingInfo().NumRings()

                print(
                    f"Reactant ring counts: {reactant_ring_counts}, Product ring count: {product_ring_count}"
                )

                # Check if any reactant has more rings than the product (forward direction)
                # OR if product has more rings than any reactant (retrosynthetic direction)
                for i, reactant_mol in enumerate(reactant_mols):
                    # Forward direction check (reactant has more rings than product)
                    if reactant_ring_counts[i] > product_ring_count:
                        reactant_smiles = reactants_smiles[i]
                        has_heterocycle = False

                        # Check for specific heterocyclic rings
                        for ring_name in heterocyclic_rings:
                            if checker.check_ring(ring_name, reactant_smiles):
                                print(f"Heterocyclic ring {ring_name} detected in reactant")
                                # Check if this ring is absent in the product
                                if not checker.check_ring(ring_name, product_smiles):
                                    print(
                                        f"Heterocyclic ring opening detected: {ring_name} ring opened"
                                    )
                                    ring_opening_detected = True
                                    has_heterocycle = True
                                    break

                        # If no specific heterocycle was identified but ring count decreased,
                        # check for heteroatoms in rings more generally
                        if not has_heterocycle:
                            # Check if the reactant has rings with heteroatoms
                            for ring in reactant_mol.GetRingInfo().AtomRings():
                                ring_has_heteroatom = False
                                for atom_idx in ring:
                                    atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                    if atom.GetSymbol() not in ["C", "H"]:
                                        ring_has_heteroatom = True
                                        print(
                                            f"Generic heterocyclic ring with {atom.GetSymbol()} detected in reactant"
                                        )
                                        break

                                if ring_has_heteroatom:
                                    # Check if this specific ring is broken in the product
                                    # This is a simplification - ideally we would use atom mapping
                                    # to track the exact atoms, but this is a reasonable approximation
                                    ring_opening_detected = True
                                    print("Generic heterocyclic ring opening detected")
                                    break

                    # Retrosynthetic direction check (product has more rings than reactant)
                    # This means in the forward direction, a ring was opened
                    elif product_ring_count > reactant_ring_counts[i]:
                        # Check if product has heterocyclic rings not in this reactant
                        for ring_name in heterocyclic_rings:
                            if checker.check_ring(
                                ring_name, product_smiles
                            ) and not checker.check_ring(ring_name, reactants_smiles[i]):
                                print(
                                    f"Retrosynthetic heterocyclic ring formation detected: {ring_name} ring formed"
                                )
                                print(f"This indicates a ring opening in the forward direction")
                                ring_opening_detected = True
                                break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Ring opening detected: {ring_opening_detected}")
    return ring_opening_detected
