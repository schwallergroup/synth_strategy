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
    This function detects if the synthetic route involves a significant cyclization in the final step.
    """
    late_cyclization_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_cyclization_detected

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if product_smiles and all(r for r in reactants_smiles):
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Count rings in product
                    if product_mol:
                        product_ring_info = product_mol.GetRingInfo()
                        product_ring_count = product_ring_info.NumRings()

                        # Count rings in reactants
                        reactant_ring_count = 0
                        reactant_mols = []
                        for reactant in reactants_smiles:
                            if reactant.strip():
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    reactant_mols.append(reactant_mol)
                                    reactant_ring_info = reactant_mol.GetRingInfo()
                                    reactant_ring_count += reactant_ring_info.NumRings()

                        # Check if product has more rings than reactants
                        if product_ring_count > reactant_ring_count:
                            # Check for known cyclization reaction types
                            rxn_smiles = node["metadata"].get("smiles", "")

                            # Check for common cyclization reactions
                            cyclization_reactions = [
                                "Paal-Knorr pyrrole synthesis",
                                "Pictet-Spengler",
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
                                "triaryl-imidazole",
                                "Fischer indole",
                                "Friedlaender chinoline",
                                "benzofuran",
                                "benzothiophene",
                                "indole",
                                "oxadiazole",
                                "Formation of NOS Heterocycles",
                                "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                                "Intramolecular amination (heterocycle formation)",
                            ]

                            is_cyclization_reaction = False
                            for rxn_type in cyclization_reactions:
                                if checker.check_reaction(rxn_type, rsmi):
                                    print(f"Detected cyclization reaction: {rxn_type}")
                                    is_cyclization_reaction = True
                                    break

                            # Check for newly formed rings
                            new_rings = []
                            for ring_name in [
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
                                "indole",
                                "quinoline",
                                "isoquinoline",
                                "benzimidazole",
                                "benzoxazole",
                                "benzothiazole",
                            ]:
                                if checker.check_ring(ring_name, product_smiles):
                                    # Check if this ring was not in any reactant
                                    ring_in_reactants = False
                                    for reactant in reactants_smiles:
                                        if checker.check_ring(ring_name, reactant):
                                            ring_in_reactants = True
                                            break

                                    if not ring_in_reactants:
                                        new_rings.append(ring_name)

                            # Check for intramolecular cyclization (one reactant forms a product with more rings)
                            intramolecular_cyclization = False
                            if (
                                len(reactants_smiles) == 1
                                and product_ring_count > reactant_ring_count
                            ):
                                intramolecular_cyclization = True
                                print("Detected intramolecular cyclization")

                            # If we have more rings in product and either a cyclization reaction or new rings formed
                            if is_cyclization_reaction or new_rings or intramolecular_cyclization:
                                late_cyclization_detected = True
                                print(
                                    f"Late-stage cyclization detected: Product has {product_ring_count} rings, reactants have {reactant_ring_count} rings"
                                )
                                if new_rings:
                                    print(f"Newly formed rings: {', '.join(new_rings)}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage cyclization strategy detected: {late_cyclization_detected}")
    return late_cyclization_detected
