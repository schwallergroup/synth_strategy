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


# Refactoring for Enumeration: Isolate the lists of chemical entities.
HETEROCYCLES_OF_INTEREST = [
    "benzofuran", "benzothiophene", "indole", "benzoxazole", "benzothiazole",
    "benzimidazole", "furan", "pyrrole", "thiophene", "oxazole", "thiazole",
    "imidazole", "pyridine", "pyrimidine", "pyrazine", "triazole", "tetrazole",
]

AMIDATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a specific heterocyclic ring is formed,
    followed by a late-stage amidation reaction. The heterocycles checked are defined
    in HETEROCYCLES_OF_INTEREST, and the amidation reactions are defined in
    AMIDATION_REACTIONS.
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
    heterocycle_formation_found = False
    late_stage_amidation = False

    # Track the depths where key reactions occur
    heterocycle_formation_depth = None
    amidation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found, late_stage_amidation
        nonlocal heterocycle_formation_depth, amidation_depth, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for heterocycle formation
                product_has_heterocycle = False
                heterocycle_formed_name = None
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_has_heterocycle = True
                        heterocycle_formed_name = heterocycle
                        break

                if product_has_heterocycle:
                    # Check if any reactant has the same heterocycle
                    reactants_have_heterocycle = False
                    for reactant in reactants_smiles:
                        for heterocycle in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(heterocycle, reactant):
                                reactants_have_heterocycle = True
                                break
                        if reactants_have_heterocycle:
                            break

                    # If heterocycle is in product but not in reactants, it was formed
                    if not reactants_have_heterocycle:
                        heterocycle_formation_found = True
                        heterocycle_formation_depth = depth
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        if heterocycle_formed_name and heterocycle_formed_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle_formed_name)
                        print(f"Heterocycle formation detected at depth {depth}")

                # Check for amidation reaction
                for reaction_type in AMIDATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Amidation reaction detected at depth {depth}: {reaction_type}")
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Track the depth of amidation (lower depth = later in synthesis)
                        if amidation_depth is None or depth < amidation_depth:
                            amidation_depth = depth
                        break

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if amidation occurs at a late stage (depth 0, 1, or 2)
    if amidation_depth is not None and amidation_depth <= 2:
        late_stage_amidation = True
        if {"type": "positional", "details": {"target": "amidation_reaction", "position_type": "depth", "operator": "<=", "value": 2, "description": "The amidation reaction must occur at a late stage in the synthesis (retrosynthetic depth of 2 or less)."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amidation_reaction", "position_type": "depth", "operator": "<=", "value": 2, "description": "The amidation reaction must occur at a late stage in the synthesis (retrosynthetic depth of 2 or less)."}})
        print(f"Late-stage amidation detected (at depth {amidation_depth})")

    # Check if heterocycle formation happens before amidation in the synthesis
    # In retrosynthetic analysis, higher depth = earlier in synthesis
    correct_sequence = (
        heterocycle_formation_depth is not None
        and amidation_depth is not None
        and heterocycle_formation_depth > amidation_depth
    )

    # The strategy is present if all conditions are met
    strategy_present = (
        heterocycle_formation_found and late_stage_amidation and correct_sequence
    )

    if heterocycle_formation_found and late_stage_amidation:
        if {"type": "co-occurrence", "details": {"targets": ["ring_formation", "amidation_reaction"], "description": "The strategy requires both the formation of a heterocycle and an amidation reaction to occur in the route."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "amidation_reaction"], "description": "The strategy requires both the formation of a heterocycle and an amidation reaction to occur in the route."}})

    if correct_sequence:
        if {"type": "sequence", "details": {"before": "ring_formation", "after": "amidation_reaction", "description": "The heterocycle formation must occur earlier in the synthesis (higher retrosynthetic depth) than the amidation reaction."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "ring_formation", "after": "amidation_reaction", "description": "The heterocycle formation must occur earlier in the synthesis (higher retrosynthetic depth) than the amidation reaction."}})

    print(f"Heterocycle formation followed by amidation strategy: {strategy_present}")
    print(
        f"Heterocycle formation depth: {heterocycle_formation_depth}, Amidation depth: {amidation_depth}"
    )

    return strategy_present, findings_json