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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a 'build-activate-couple' strategy where a complex core is built,
    then activated (via functional group transformations) for final coupling.
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

    # Track functional group transformations
    fg_transformations = []

    # Track if we found key steps
    found_coupling = False
    found_activation = False
    found_core_building = False

    def dfs_traverse(node, depth=0):
        nonlocal found_coupling, found_activation, found_core_building
        nonlocal fg_transformations, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")

            # Check for coupling reactions (final coupling step)
            if depth <= 3:  # Late stage - expanded depth range
                # Check for various coupling reactions
                coupling_reactions = [
                    "Schotten-Baumann to ester",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with secondary amine to amide",
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Stille reaction_aryl",
                    "Negishi coupling",
                    "Heck_terminal_vinyl",
                    "Buchwald-Hartwig",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Sonogashira acetylene_aryl halide",
                    "Sonogashira alkyne_aryl halide",
                ]

                for rxn_name in coupling_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        found_coupling = True
                        fg_transformations.append(("coupling", depth))
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        print(f"Found coupling reaction ({rxn_name}) at depth {depth}")
                        break

            # Check for activation (functional group transformations)
            if depth <= 6:  # Mid-stage - expanded depth range
                # Check for various activation reactions
                activation_reactions = [
                    "Acyl chlorides from alcohols",
                    "Acyl bromides from alcohols",
                    "Acyl iodides from alcohols",
                    "Carboxylic acid to amide conversion",
                    "Esterification of Carboxylic Acids",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "COOH ethyl deprotection",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to triflate conversion",
                    "Preparation of boronic acids",
                    "Preparation of boronic esters",
                    "Aromatic halogenation",
                    "Aromatic bromination",
                    "Aromatic chlorination",
                    "Aromatic iodination",
                    "Aromatic fluorination",
                    "Oxidation of aldehydes to carboxylic acids",
                    "Oxidation of alcohol to carboxylic acid",
                    "Oxidation of ketone to carboxylic acid",
                ]

                for rxn_name in activation_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        found_activation = True
                        fg_transformations.append(("activation", depth))
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        print(f"Found activation reaction ({rxn_name}) at depth {depth}")
                        break

            # Check for core building (cyclization, ring formation)
            if depth >= 2:  # Early stage - relaxed depth requirement
                # Check for specific cyclization reactions - expanded list
                cyclization_reactions = [
                    "Formation of NOS Heterocycles",
                    "Intramolecular transesterification/Lactone formation",
                    "Diels-Alder",
                    "Paal-Knorr pyrrole synthesis",
                    "Intramolecular amination (heterocycle formation)",
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                    "Pictet-Spengler",
                    "Fischer indole",
                    "Friedlaender chinoline",
                    "benzimidazole_derivatives_carboxylic-acid/ester",
                    "benzimidazole_derivatives_aldehyde",
                    "benzothiazole",
                    "benzoxazole_arom-aldehyde",
                    "benzoxazole_carboxylic-acid",
                    "thiazole",
                    "tetrazole_terminal",
                    "Huisgen_Cu-catalyzed_1,4-subst",
                    "Huisgen_Ru-catalyzed_1,5_subst",
                    "1,2,4-triazole_acetohydrazide",
                    "1,2,4-triazole_carboxylic-acid/ester",
                    "pyrazole",
                    "Paal-Knorr pyrrole",
                    "benzofuran",
                    "benzothiophene",
                    "indole",
                    "oxadiazole",
                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                    "Huisgen 1,3 dipolar cycloaddition",
                    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                ]

                for rxn_name in cyclization_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        found_core_building = True
                        fg_transformations.append(("cyclization", depth))
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        print(f"Found cyclization ({rxn_name}) at depth {depth}")
                        break

                # If no specific cyclization reaction is found, check for ring count changes
                if not found_core_building:
                    reactant_rings = 0
                    product_rings = 0

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                reactant_rings += mol.GetRingInfo().NumRings()
                        except:
                            continue

                    try:
                        product_mol = Chem.MolFromSmiles(product_part)
                        if product_mol:
                            product_rings = product_mol.GetRingInfo().NumRings()
                    except:
                        pass

                    # If product has more rings than reactants, it's a cyclization
                    if product_rings > reactant_rings:
                        found_core_building = True
                        fg_transformations.append(("cyclization", depth))
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        print(f"Found cyclization (ring count increased) at depth {depth}")

                    # Also check for common ring structures in the product
                    common_rings = [
                        "pyridine",
                        "pyrrole",
                        "furan",
                        "thiophene",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "indole",
                        "benzofuran",
                        "quinoline",
                        "isoquinoline",
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                        "triazole",
                        "tetrazole",
                        "piperidine",
                        "morpholine",
                        "thiomorpholine",
                        "pyrrolidine",
                        "piperazine",
                    ]

                    for ring_name in common_rings:
                        if checker.check_ring(ring_name, product_part):
                            # Check if this ring wasn't in the reactants
                            ring_in_reactants = False
                            for reactant in reactants:
                                if checker.check_ring(ring_name, reactant):
                                    ring_in_reactants = True
                                    break

                            if not ring_in_reactants:
                                found_core_building = True
                                fg_transformations.append(("cyclization", depth))
                                if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                                print(f"Found cyclization (formed {ring_name}) at depth {depth}")
                                break

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Found steps - Coupling: {found_coupling}, Activation: {found_activation}, Core building: {found_core_building}"
    )
    print(f"Transformations: {fg_transformations}")

    result = False
    # More flexible detection logic
    # If we found core building, that's the minimum requirement
    if found_core_building:
        # Add the 'core_building' co-occurrence constraint if found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "core_building"
                ]
            }
        })

        # Sort by depth (ascending)
        fg_transformations.sort(key=lambda x: x[1])

        # Find the minimum depth for each transformation type
        min_cyclization_depth = float("inf")
        min_activation_depth = float("inf")
        min_coupling_depth = float("inf")

        for transform_type, depth in fg_transformations:
            if transform_type in ["cyclization", "core_building"] and depth < min_cyclization_depth:
                min_cyclization_depth = depth
            elif (
                transform_type in ["activation", "activation_fg_change"]
                and depth < min_activation_depth
            ):
                min_activation_depth = depth
            elif transform_type == "coupling" and depth < min_coupling_depth:
                min_coupling_depth = depth

        print(
            f"Min depths - Core building: {min_cyclization_depth}, Activation: {min_activation_depth}, Coupling: {min_coupling_depth}"
        )

        # Core building alone is sufficient for a simplified build-activate-couple strategy
        if min_cyclization_depth != float("inf"):
            # If we have all three steps in the correct order
            if (
                min_cyclization_depth < min_activation_depth < min_coupling_depth
                and min_activation_depth != float("inf")
                and min_coupling_depth != float("inf")
            ):
                print("Detected complete build-activate-couple strategy")
                result = True
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "ordered_events": [
                            "core_building",
                            "activation",
                            "coupling"
                        ]
                    }
                })

            # If we have building and activation in the correct order
            elif min_cyclization_depth < min_activation_depth and min_activation_depth != float(
                "inf"
            ):
                print("Detected build-activate strategy (coupling step not found or implicit)")
                result = True
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "ordered_events": [
                            "core_building",
                            "activation"
                        ]
                    }
                })

            # If we have building and coupling in the correct order
            elif min_cyclization_depth < min_coupling_depth and min_coupling_depth != float("inf"):
                print("Detected build-couple strategy (activation step not found or implicit)")
                result = True
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "ordered_events": [
                            "core_building",
                            "coupling"
                        ]
                    }
                })

            # If we only have core building, that's still a simplified version of the strategy
            else:
                print(
                    "Detected simplified build strategy (activation and coupling steps not found or implicit)"
                )
                result = True

    # Add positional constraints if the corresponding steps were found
    if found_coupling:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "coupling",
                "position": "late_stage",
                "max_depth": 3
            }
        })
    if found_activation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "activation",
                "position": "mid_stage",
                "max_depth": 6
            }
        })
    if found_core_building:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "core_building",
                "position": "early_stage",
                "min_depth": 2
            }
        })

    return result, findings_json
