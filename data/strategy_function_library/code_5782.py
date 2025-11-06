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


HETEROCYCLES_OF_INTEREST = [
    "benzoxazole",
    "benzimidazole",
    "benzothiazole",
    "oxazole",
    "imidazole",
    "thiazole",
    "pyrrole",
    "pyrazole",
    "triazole",
    "tetrazole",
    "indole",
    "quinoline",
    "isoquinoline",
    "purine",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acylation of primary amines",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles",
    "Benzoxazole formation from aldehyde",
    "Benzoxazole formation from acyl halide",
    "Benzoxazole formation from ester/carboxylic acid",
    "Benzoxazole formation (intramolecular)",
    "Benzimidazole formation from aldehyde",
    "Benzimidazole formation from acyl halide",
    "Benzimidazole formation from ester/carboxylic acid",
    "Benzothiazole formation from aldehyde",
    "Benzothiazole formation from acyl halide",
    "Benzothiazole formation from ester/carboxylic acid",
    "Paal-Knorr pyrrole synthesis",
    "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzimidazole_derivatives_aldehyde}",
    "{benzothiazole}",
    "{benzoxazole_arom-aldehyde}",
    "{benzoxazole_carboxylic-acid}",
    "{thiazole}",
    "{tetrazole_terminal}",
    "{1,2,4-triazole_acetohydrazide}",
    "{1,2,4-triazole_carboxylic-acid/ester}",
    "{pyrazole}",
    "{Paal-Knorr pyrrole}",
    "{Fischer indole}",
    "{indole}",
    "{oxadiazole}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a specific three-step synthetic sequence: 1) Reduction of a nitro group to an amine, followed by 2) an amide formation reaction, followed by 3) the formation of a heterocycle. The specific amide formation reactions, heterocycle formation reactions, and target heterocycles are defined in the module-level constants AMIDE_FORMATION_REACTIONS, HETEROCYCLE_FORMATION_REACTIONS, and HETEROCYCLES_OF_INTEREST, respectively. The function identifies the reaction steps and verifies they occur in the correct retrosynthetic order (heterocycle formation at the lowest depth).
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

    # Initialize sequence tracking
    nitro_reduction_depth = None
    amide_formation_depth = None
    heterocycle_formation_depth = None

    # Track molecules with specific functional groups
    nitro_molecules = {}  # {mol_smiles: depth}
    amine_molecules = {}  # {mol_smiles: depth}
    amide_molecules = {}  # {mol_smiles: depth}
    heterocycle_molecules = {}  # {mol_smiles: depth}

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_depth, amide_formation_depth, heterocycle_formation_depth, findings_json

        # Process molecule nodes to track functional groups
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for nitro groups
            if checker.check_fg("Nitro group", mol_smiles):
                nitro_molecules[mol_smiles] = depth
                if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

            # Check for amines
            if checker.check_fg("Primary amine", mol_smiles):
                amine_molecules[mol_smiles] = depth
                if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
            if checker.check_fg("Secondary amine", mol_smiles):
                amine_molecules[mol_smiles] = depth
                if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
            if checker.check_fg("Tertiary amine", mol_smiles):
                amine_molecules[mol_smiles] = depth
                if "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")

            # Check for amides
            if checker.check_fg("Primary amide", mol_smiles):
                amide_molecules[mol_smiles] = depth
                if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
            if checker.check_fg("Secondary amide", mol_smiles):
                amide_molecules[mol_smiles] = depth
                if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
            if checker.check_fg("Tertiary amide", mol_smiles):
                amide_molecules[mol_smiles] = depth
                if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

            # Check for heterocycles
            for ring in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring, mol_smiles):
                    heterocycle_molecules[mol_smiles] = depth
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    break

        # Process reaction nodes to identify transformations
        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    has_nitro_reactant = any(
                        checker.check_fg("Nitro group", reactant) for reactant in reactants
                    )
                    has_amine_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )

                    if has_nitro_reactant and has_amine_product:
                        nitro_reduction_depth = depth
                        if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                # Check for amide formation
                for rxn in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        has_amine_reactant = any(
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            for reactant in reactants
                        )
                        has_amide_product = (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        )

                        if has_amine_reactant and has_amide_product:
                            amide_formation_depth = depth
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                # Check for heterocycle formation
                # First check specific heterocycle formation reactions
                for rxn in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        has_amide_reactant = any(
                            checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                            for reactant in reactants
                        )
                        has_heterocycle_product = any(
                            checker.check_ring(ring, product) for ring in HETEROCYCLES_OF_INTEREST
                        )

                        if has_amide_reactant and has_heterocycle_product:
                            heterocycle_formation_depth = depth
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                # If no specific reaction was found, check for general heterocycle formation
                if heterocycle_formation_depth is None:
                    has_heterocycle_product = any(
                        checker.check_ring(ring, product) for ring in HETEROCYCLES_OF_INTEREST
                    )
                    reactants_have_heterocycle = any(
                        any(checker.check_ring(ring, reactant) for ring in HETEROCYCLES_OF_INTEREST)
                        for reactant in reactants
                    )
                    has_amide_reactant = any(
                        checker.check_fg("Primary amide", reactant)
                        or checker.check_fg("Secondary amide", reactant)
                        or checker.check_fg("Tertiary amide", reactant)
                        for reactant in reactants
                    )

                    if (
                        has_heterocycle_product
                        and not reactants_have_heterocycle
                        and has_amide_reactant
                    ):
                        heterocycle_formation_depth = depth
                        # Add a generic 'ring_formation' if not already added by a specific reaction
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children nodes
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = False
    # Check if all steps were found and in the correct order
    if (
        nitro_reduction_depth is not None
        and amide_formation_depth is not None
        and heterocycle_formation_depth is not None
    ):
        correct_sequence = (
            nitro_reduction_depth > amide_formation_depth > heterocycle_formation_depth
        )
        if correct_sequence:
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "heterocycle_formation",
                        "amide_formation",
                        "nitro_reduction"
                    ],
                    "direction": "retrosynthetic"
                }
            })

    # Check if we can infer the sequence from molecule tracking
    elif len(heterocycle_molecules) > 0 and len(amide_molecules) > 0 and len(nitro_molecules) > 0:
        min_heterocycle_depth = min(heterocycle_molecules.values())
        min_amide_depth = min(amide_molecules.values())
        min_nitro_depth = min(nitro_molecules.values())

        correct_sequence = min_nitro_depth > min_amide_depth > min_heterocycle_depth
        if correct_sequence:
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "heterocycle_formation",
                        "amide_formation",
                        "nitro_reduction"
                    ],
                    "direction": "retrosynthetic"
                }
            })

    return result, findings_json
