from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

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


HYDRAZINE_DERIVED_HETEROCYCLES = ['pyrazole', 'triazole', 'tetrazole', 'indazole']

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy where hydrazine is incorporated via
    nucleophilic aromatic substitution and later used for heterocycle formation.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track the key steps in the strategy
    hydrazine_incorporation_reactions = []
    heterocycle_formation_reactions = []

    def dfs_traverse(node, depth=0, path=[]):
        nonlocal hydrazine_incorporation_reactions, heterocycle_formation_reactions, findings_json
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydrazine incorporation via nucleophilic substitution
                hydrazine_reactant = False

                for reactant in reactants:
                    if checker.check_fg("Hydrazine", reactant):
                        hydrazine_reactant = True
                        if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")
                    if checker.check_fg("Acylhydrazine", reactant):
                        hydrazine_reactant = True
                        if "Acylhydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acylhydrazine")
                    if checker.check_fg("Hydrazone", reactant):
                        hydrazine_reactant = True
                        if "Hydrazone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Hydrazone")

                # Check for hydrazine incorporation reactions
                if hydrazine_reactant:
                    # Check for nucleophilic substitution reactions that incorporate hydrazine
                    nuc_sub_reactions = [
                        "heteroaromatic_nuc_sub",
                        "nucl_sub_aromatic_ortho_nitro",
                        "nucl_sub_aromatic_para_nitro",
                        "N-arylation",
                        "Buchwald-Hartwig",
                        "N-arylation_heterocycles"
                    ]
                    for rxn_name in nuc_sub_reactions:
                        if checker.check_reaction(rxn_name, rsmi):
                            print(
                                f"Hydrazine incorporation via nucleophilic aromatic substitution detected at depth {depth}"
                            )
                            hydrazine_incorporation_reactions.append((depth, rsmi))
                            if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                            break # Only need one to match

                # Check if hydrazine is used for heterocycle formation
                hydrazine_in_reactants = any(
                    checker.check_fg("Hydrazine", r)
                    or checker.check_fg("Acylhydrazine", r)
                    or checker.check_fg("Hydrazone", r)
                    for r in reactants
                )
                heterocycle_in_product = False
                for ring in HYDRAZINE_DERIVED_HETEROCYCLES:
                    if checker.check_ring(ring, product):
                        heterocycle_in_product = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

                if hydrazine_in_reactants and heterocycle_in_product:
                    # Check for specific heterocycle formation reactions if possible
                    heterocycle_formation_rxns = [
                        "pyrazole",
                        "[3+2]-cycloaddition of hydrazone and alkyne",
                        "[3+2]-cycloaddition of hydrazone and alkene",
                        "Huisgen"
                    ]
                    found_specific_heterocycle_rxn = False
                    for rxn_name in heterocycle_formation_rxns:
                        if checker.check_reaction(rxn_name, rsmi):
                            print(f"Hydrazine used for heterocycle formation at depth {depth}")
                            heterocycle_formation_reactions.append((depth, rsmi))
                            if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                            found_specific_heterocycle_rxn = True
                            break

                    if not found_specific_heterocycle_rxn:
                        # If no specific reaction check passes, still consider it if we have heterocycle in product
                        print(f"Potential hydrazine heterocycle formation at depth {depth}")
                        heterocycle_formation_reactions.append((depth, rsmi))
                        # Add a generic 'ring_formation' if not already present and no specific reaction was found
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # This means it's a 'chemical' node or similar
                new_depth = depth + 1
            dfs_traverse(child, new_depth, path + [node])

    # Start traversal from the root
    dfs_traverse(route)

    result = False

    # Check if we found both steps in the correct sequence
    if not hydrazine_incorporation_reactions and not heterocycle_formation_reactions:
        print("No hydrazine incorporation or heterocycle formation reactions found")
        return result, findings_json

    # If we only found heterocycle formation but not incorporation, check if hydrazine might be a starting material
    if heterocycle_formation_reactions and not hydrazine_incorporation_reactions:
        print(
            "Found heterocycle formation but no explicit hydrazine incorporation - checking for hydrazine as starting material"
        )

        # Check if hydrazine or derivatives are in-stock materials
        def check_starting_materials(node):
            nonlocal findings_json
            if node["type"] == "mol" and node.get("in_stock", False):
                if checker.check_fg("Hydrazine", node["smiles"]):
                    if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")
                    return True
                if checker.check_fg("Acylhydrazine", node["smiles"]):
                    if "Acylhydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Acylhydrazine")
                    return True
                if checker.check_fg("Hydrazone", node["smiles"]):
                    if "Hydrazone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Hydrazone")
                    return True
            for child in node.get("children", []):
                if check_starting_materials(child):
                    return True
            return False

        if check_starting_materials(route):
            print("Hydrazine or derivative found as starting material")
            # Assign a very high depth to the incorporation (earlier in synthesis)
            hydrazine_incorporation_reactions.append((float("inf"), "starting_material"))
            # Record the structural constraint for starting material
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "description": "A hydrazine moiety is first introduced into a molecule (or present as a starting material) and is subsequently used to construct a heterocycle.",
                    "event_before_group": {
                        "description": "The incorporation of a hydrazine derivative, satisfied by any of the following reactions or conditions.",
                        "reactions": [],
                        "conditions": [
                            "hydrazine_derivative_as_starting_material"
                        ]
                    },
                    "event_after_group": {
                        "description": "The formation of a hydrazine-derived heterocycle, satisfied by any of the following reactions.",
                        "reactions": [
                            "pyrazole",
                            "[3+2]-cycloaddition of hydrazone and alkyne",
                            "[3+2]-cycloaddition of hydrazone and alkene",
                            "Huisgen",
                            "ring_formation"
                        ],
                        "conditions": []
                    }
                }
            })

    # Find the minimum depth for each reaction type
    min_incorporation_depth = (
        min([d for d, _ in hydrazine_incorporation_reactions])
        if hydrazine_incorporation_reactions
        else float("inf")
    )
    min_heterocycle_depth = (
        min([d for d, _ in heterocycle_formation_reactions])
        if heterocycle_formation_reactions
        else float("inf")
    )

    # The incorporation should happen before heterocycle formation in the synthesis direction
    # In retrosynthesis traversal, this means incorporation depth > heterocycle depth
    print(
        f"Min incorporation depth: {min_incorporation_depth}, Min heterocycle depth: {min_heterocycle_depth}"
    )

    # If we have both steps, check the sequence
    if hydrazine_incorporation_reactions and heterocycle_formation_reactions:
        if min_incorporation_depth > min_heterocycle_depth:
            result = True
            # Record the structural constraint for sequence
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "description": "A hydrazine moiety is first introduced into a molecule (or present as a starting material) and is subsequently used to construct a heterocycle.",
                    "event_before_group": {
                        "description": "The incorporation of a hydrazine derivative, satisfied by any of the following reactions or conditions.",
                        "reactions": [
                            "heteroaromatic_nuc_sub",
                            "nucl_sub_aromatic_ortho_nitro",
                            "nucl_sub_aromatic_para_nitro",
                            "N-arylation",
                            "Buchwald-Hartwig",
                            "N-arylation_heterocycles"
                        ],
                        "conditions": []
                    },
                    "event_after_group": {
                        "description": "The formation of a hydrazine-derived heterocycle, satisfied by any of the following reactions.",
                        "reactions": [
                            "pyrazole",
                            "[3+2]-cycloaddition of hydrazone and alkyne",
                            "[3+2]-cycloaddition of hydrazone and alkene",
                            "Huisgen",
                            "ring_formation"
                        ],
                        "conditions": []
                    }
                }
            })

    # If we only have heterocycle formation and hydrazine is a starting material
    if heterocycle_formation_reactions and min_incorporation_depth == float("inf"):
        print("Strategy detected with hydrazine as starting material")
        result = True

    return result, findings_json