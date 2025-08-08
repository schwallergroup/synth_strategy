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
    This function detects late-stage formylation on a heterocyclic scaffold.
    Late-stage means the formylation occurs in the final synthetic step (depth 0).
    """
    print("Starting late_stage_formylation_strategy analysis")
    formylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal formylation_detected

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                try:
                    print(f"Analyzing late-stage reaction: {rsmi}")

                    # Check if product has an aldehyde group
                    has_aldehyde_in_product = checker.check_fg("Aldehyde", product)
                    print(f"Product has aldehyde: {has_aldehyde_in_product}")

                    if not has_aldehyde_in_product:
                        print("No aldehyde in product, skipping")
                        return  # No aldehyde in product, can't be formylation

                    # Check if product contains heterocyclic rings
                    heterocyclic_rings = [
                        "furan",
                        "pyrrole",
                        "pyridine",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrimidine",
                        "pyrazine",
                        "thiophene",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                        "triazole",
                        "tetrazole",
                        "oxadiazole",
                        "thiadiazole",
                        "isoxazole",
                        "isothiazole",
                    ]

                    has_heterocycle = False
                    for ring in heterocyclic_rings:
                        if checker.check_ring(ring, product):
                            has_heterocycle = True
                            print(f"Found heterocyclic ring: {ring}")
                            break

                    if not has_heterocycle:
                        print("No heterocycle in product, skipping")
                        return  # No heterocycle in product, not a match

                    # Check if any reactant has aldehyde group and count aldehydes
                    has_aldehyde_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Aldehyde", reactant):
                            has_aldehyde_in_reactants = True
                            print(f"Found aldehyde in reactant: {reactant}")
                            break

                    # Check for formylation reactions directly
                    formylation_reactions = [
                        "Friedel-Crafts acylation",
                        "Aromatic hydroxylation",
                        "Bouveault aldehyde synthesis",
                        "Oxidation of alkene to aldehyde",
                        "Oxidation of primary alcohol to aldehyde",
                        "Duff reaction",
                        "Reimer-Tiemann reaction",
                        "Vilsmeier-Haack reaction",
                    ]

                    is_formylation_reaction = False
                    for rxn_type in formylation_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_formylation_reaction = True
                            print(f"Detected formylation reaction: {rxn_type}")
                            break

                    # Check for common formylating agents
                    formylating_agents = [
                        "CN(C)C=O",  # DMF
                        "O=CH",  # Formaldehyde
                        "O=CHO",  # Formic acid
                        "O=COC=O",  # Formic anhydride
                        "ClC=O",  # Formyl chloride
                        "O=CN",  # Formamide
                    ]

                    is_formylation = False
                    for reactant in reactants:
                        # Check for formaldehyde specifically
                        if checker.check_fg("Formaldehyde", reactant):
                            print(f"Found formaldehyde in reactant: {reactant}")
                            is_formylation = True
                            break

                        # Check for other formylating agents
                        for agent in formylating_agents:
                            try:
                                agent_mol = Chem.MolFromSmiles(agent)
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if (
                                    agent_mol
                                    and reactant_mol
                                    and reactant_mol.HasSubstructMatch(agent_mol)
                                ):
                                    is_formylation = True
                                    print(f"Found formylating agent: {agent} in {reactant}")
                                    break
                            except Exception as e:
                                print(f"Error checking formylating agent: {e}")
                                continue

                        if is_formylation:
                            break

                    # Check for Vilsmeier-Haack formylation specifically
                    has_vilsmeier = False
                    has_dmf = False

                    for reactant in reactants:
                        # Check for DMF
                        if (
                            "CN(C)C=O" in reactant
                            or "OC=N(C)C" in reactant
                            or "CN(C)CHO" in reactant
                        ):
                            has_dmf = True
                            print(f"Found DMF in reactant: {reactant}")

                        # Check for phosphorus reagents
                        if (
                            checker.check_fg("Phosphate ester", reactant)
                            or "POCl3" in reactant
                            or "PCl3" in reactant
                            or "PCl5" in reactant
                        ):
                            has_vilsmeier = True
                            print(f"Found phosphorus reagent in reactant: {reactant}")

                    if has_dmf and has_vilsmeier:
                        is_formylation = True
                        print("Found Vilsmeier-Haack formylation conditions")

                    # Formylation detected if:
                    # 1. Product has aldehyde group
                    # 2. Reactants don't have aldehyde group (it's being added) OR it's a known formylation reaction
                    # 3. Product has heterocyclic scaffold
                    # 4. Formylating agent or formylation reaction is present
                    if (
                        has_aldehyde_in_product
                        and has_heterocycle
                        and (not has_aldehyde_in_reactants or is_formylation_reaction)
                        and (is_formylation or is_formylation_reaction)
                    ):
                        print(f"Late-stage formylation detected on heterocyclic scaffold")
                        formylation_detected = True

                except Exception as e:
                    print(f"Error processing SMILES in formylation detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: formylation_detected = {formylation_detected}")
    return formylation_detected
