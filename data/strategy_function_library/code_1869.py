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


HETEROCYCLE_RINGS_OF_INTEREST = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazine", "quinoline",
    "isoquinoline", "quinazoline", "quinoxaline", "phthalazine", "furan",
    "thiophene", "pyrrole", "oxazole", "thiazole", "imidazole", "triazole",
    "tetrazole", "piperidine", "piperazine", "morpholine", "indole",
    "benzimidazole", "benzoxazole", "benzothiazole", "purine",
]

AMINE_SUBSTITUTION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "reductive amination",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
    "Goldberg coupling",
    "Ullmann-Goldberg Substitution amine",
    "aza-Michael addition primary",
    "aza-Michael addition secondary",
    "aza-Michael addition aromatic",
    "Ring opening of epoxide with amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy where a haloalkoxy linker is installed on an aromatic core,
    followed by scaffold-building reactions, and finally a late-stage amine substitution.
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

    # Track if we found the key features
    found_linker_installation = False
    found_heterocycle_formation = False
    found_late_stage_amine_substitution = False

    # Track the depth at which each feature was found
    linker_installation_depth = -1
    heterocycle_formation_depth = -1
    amine_substitution_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_linker_installation, found_heterocycle_formation, found_late_stage_amine_substitution
        nonlocal linker_installation_depth, heterocycle_formation_depth, amine_substitution_depth
        nonlocal findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for linker installation (broader definition)
                # Linker installation can be via ether formation, ester formation, or amide formation
                linker_check_passed = False
                if (
                    (
                        any(checker.check_fg("Phenol", r) for r in reactants_smiles)
                        and (
                            any(
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Aromatic halide", r)
                                or checker.check_fg("Triflate", r)
                                or checker.check_fg("Mesylate", r)
                                or checker.check_fg("Tosylate", r)
                                for r in reactants_smiles
                            )
                        )
                        and (
                            checker.check_reaction("Williamson Ether Synthesis", rsmi)
                            or checker.check_reaction("Mitsunobu_phenole", rsmi)
                            or checker.check_reaction("Williamson ether", rsmi)
                        )
                    )
                ):
                    linker_check_passed = True
                    for r in reactants_smiles:
                        if checker.check_fg("Phenol", r):
                            if "Phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                        if checker.check_fg("Primary halide", r):
                            if "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
                        if checker.check_fg("Secondary halide", r):
                            if "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
                        if checker.check_fg("Aromatic halide", r):
                            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        if checker.check_fg("Triflate", r):
                            if "Triflate" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Triflate")
                        if checker.check_fg("Mesylate", r):
                            if "Mesylate" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
                        if checker.check_fg("Tosylate", r):
                            if "Tosylate" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tosylate")
                    if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                        if "Williamson Ether Synthesis" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")
                    if checker.check_reaction("Mitsunobu_phenole", rsmi):
                        if "Mitsunobu_phenole" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu_phenole")
                    if checker.check_reaction("Williamson ether", rsmi):
                        if "Williamson ether" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Williamson ether")
                
                if (
                    any(checker.check_fg("Phenol", r) for r in reactants_smiles)
                    and checker.check_fg("Ether", product_smiles)
                    and not all(checker.check_fg("Ether", r) for r in reactants_smiles)
                ):
                    linker_check_passed = True
                    if "Phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                    if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ether")

                if (
                    any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                    and any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants_smiles
                    )
                    and checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                ):
                    linker_check_passed = True
                    for r in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", r):
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if checker.check_fg("Primary alcohol", r):
                            if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                        if checker.check_fg("Secondary alcohol", r):
                            if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                        if checker.check_fg("Tertiary alcohol", r):
                            if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                    if "Esterification of Carboxylic Acids" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")

                if (
                    any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                    and any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants_smiles
                    )
                    and (
                        checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        )
                    )
                ):
                    linker_check_passed = True
                    for r in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", r):
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if checker.check_fg("Primary amine", r):
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", r):
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                    if checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi):
                        if "Carboxylic acid with primary amine to amide" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Carboxylic acid with primary amine to amide")
                    if checker.check_reaction("Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi):
                        if "Acylation of Nitrogen Nucleophiles by Carboxylic Acids" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Acylation of Nitrogen Nucleophiles by Carboxylic Acids")

                if linker_check_passed:
                    found_linker_installation = True
                    linker_installation_depth = depth

                # Check for heterocycle formation
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

                if all([product_mol] + reactant_mols):  # Ensure all molecules parsed correctly
                    product_ring_count = len(Chem.GetSSSR(product_mol))
                    reactant_ring_count = sum(len(Chem.GetSSSR(r)) for r in reactant_mols)

                    if product_ring_count > reactant_ring_count:
                        # Check if any of the new rings are heterocycles
                        for ring in HETEROCYCLE_RINGS_OF_INTEREST:
                            if checker.check_ring(ring, product_smiles) and not any(
                                checker.check_ring(ring, r) for r in reactants_smiles
                            ):
                                found_heterocycle_formation = True
                                heterocycle_formation_depth = depth
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                break

                # Check for late-stage amine substitution
                amine_substitution_check_passed = False
                for rxn in AMINE_SUBSTITUTION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        amine_substitution_check_passed = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                if amine_substitution_check_passed:
                    found_late_stage_amine_substitution = True
                    amine_substitution_depth = depth

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = False
    # Check for the strategy with correct ordering
    if found_heterocycle_formation:
        if found_linker_installation and found_late_stage_amine_substitution:
            if (
                linker_installation_depth > heterocycle_formation_depth
                and heterocycle_formation_depth > amine_substitution_depth
            ):
                result = True
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["linker_installation", "heterocycle_formation"], "description": "A linker installation step must occur earlier in the synthesis (higher depth) than a heterocycle formation step."}})
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["heterocycle_formation", "late_stage_amine_substitution"], "description": "A heterocycle formation step must occur earlier in the synthesis (higher depth) than a late-stage amine substitution."}})
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["linker_installation", "late_stage_amine_substitution"], "description": "A linker installation step must occur earlier in the synthesis (higher depth) than a late-stage amine substitution."}})
        elif found_late_stage_amine_substitution:
            if heterocycle_formation_depth > amine_substitution_depth:
                result = True
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["heterocycle_formation", "late_stage_amine_substitution"], "description": "A heterocycle formation step must occur earlier in the synthesis (higher depth) than a late-stage amine substitution."}})
        elif found_linker_installation:
            if linker_installation_depth > heterocycle_formation_depth:
                result = True
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["linker_installation", "heterocycle_formation"], "description": "A linker installation step must occur earlier in the synthesis (higher depth) than a heterocycle formation step."}})

    if (
        found_linker_installation
        and found_late_stage_amine_substitution
        and not found_heterocycle_formation
    ):
        if linker_installation_depth > amine_substitution_depth:
            result = True
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["linker_installation", "late_stage_amine_substitution"], "description": "A linker installation step must occur earlier in the synthesis (higher depth) than a late-stage amine substitution."}})
            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "heterocycle_formation", "description": "One valid synthetic path requires that no heterocycle formation occurs throughout the route."}})

    return result, findings_json