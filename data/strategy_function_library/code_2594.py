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


HETEROAROMATIC_RINGS_FOR_SNAR = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "imidazole", "oxazole", "thiazole",
]

SNAR_LIKE_REACTION_TYPES = [
    "heteroaromatic_nuc_sub", "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro", "N-arylation_heterocycles", "Buchwald-Hartwig",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects sequential SNAr reactions on specific halogenated heteroaromatics, followed by late-stage hydrazide formation.
    The checked heteroaromatics and reaction types are defined in HETEROAROMATIC_RINGS_FOR_SNAR and SNAR_LIKE_REACTION_TYPES, respectively.
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

    snar_reactions = []
    hydrazide_formation_reaction = None
    late_stage_depth = float("inf")

    def is_heteroaromatic_with_halide(smiles):
        """Check if molecule is a heteroaromatic with a halide"""
        is_aromatic_halide = checker.check_fg("Aromatic halide", smiles)
        if is_aromatic_halide:
            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

        found_ring = False
        for ring in HETEROAROMATIC_RINGS_FOR_SNAR:
            if checker.check_ring(ring, smiles):
                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                found_ring = True

        return is_aromatic_halide and found_ring

    def dfs_traverse(node, depth=0):
        nonlocal snar_reactions, hydrazide_formation_reaction, late_stage_depth, findings_json

        # Track the minimum depth to identify late-stage reactions
        if node["type"] == "reaction" and depth < late_stage_depth:
            late_stage_depth = depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr reaction on heteroaromatics
            for reactant in reactants:
                if is_heteroaromatic_with_halide(reactant):
                    # Check if this is an SNAr reaction
                    for rxn_type in SNAR_LIKE_REACTION_TYPES:
                        if checker.check_reaction(rxn_type, rsmi):
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            print(f"Found heteroaromatic SNAr reaction: {rsmi}")
                            snar_reactions.append((rsmi, depth, product))
                            break
                    break

            # Check for hydrazide formation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and checker.check_fg("Acylhydrazine", product):
                if "Acylhydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acylhydrazine")

                # Verify acylhydrazine is being formed, not just present
                is_formed = not all(checker.check_fg("Acylhydrazine", r) for r in reactants)
                if is_formed:
                    if "hydrazide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("hydrazide_formation")
                    print(f"Found hydrazide formation: {rsmi}")
                    # Store the reaction and its depth
                    if not hydrazide_formation_reaction or depth < hydrazide_formation_reaction[1]:
                        hydrazide_formation_reaction = (rsmi, depth)
                    # Add negation constraint if applicable
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Acylhydrazine", "scope": "reactants", "description": "To identify a true 'hydrazide_formation' reaction, the Acylhydrazine group must not be present in all reactants of the step that forms it."}})

        # Traverse children with increased depth
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    def are_reactions_sequential(reactions):
        """
        Check if the SNAr reactions are sequential (building on the same molecule)
        Uses common heteroaromatic cores to track molecules through the synthesis
        """
        if len(reactions) < 2:
            return False

        # Sort reactions by depth (from early to late stage)
        sorted_reactions = sorted(reactions, key=lambda x: x[1], reverse=True)

        # Track if we found at least one sequential pair
        found_sequential = False

        for i in range(len(sorted_reactions) - 1):
            current_rsmi, _, current_product = sorted_reactions[i]
            current_product_mol = Chem.MolFromSmiles(current_product)

            if not current_product_mol:
                continue

            for j in range(i + 1, len(sorted_reactions)):
                next_rsmi, _, _ = sorted_reactions[j]
                next_reactants = next_rsmi.split(">")[0].split(".")

                for reactant in next_reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Check if the reactant contains a significant part of the previous product
                    # by checking for common heteroaromatic core
                    for ring in HETEROAROMATIC_RINGS_FOR_SNAR:
                        if checker.check_ring(ring, current_product) and checker.check_ring(
                            ring, reactant
                        ):
                            print(
                                f"Sequential reactions found: product {current_product} -> reactant {reactant}"
                            )
                            found_sequential = True
                            break

                    if found_sequential:
                        break

                if found_sequential:
                    break

            if found_sequential:
                break

        if found_sequential:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"description": "The SNAr reactions must be sequential, where the product of one (identified by its heteroaromatic core) serves as a reactant for a subsequent one.", "ordered_events": ["SNAr reaction on specified heteroaromatic", "SNAr reaction on specified heteroaromatic"], "on_same_scaffold": True}})

        return found_sequential

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple SNAr reactions
    has_multiple_snar = len(snar_reactions) >= 2
    if has_multiple_snar:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "SNAr reaction on specified heteroaromatic", "operator": ">=", "value": 2, "description": "Requires at least two SNAr-like reactions on a reactant containing both an aromatic halide and a specified heteroaromatic ring."}})

    # Check if the SNAr reactions are sequential
    are_sequential = are_reactions_sequential(snar_reactions) if has_multiple_snar else False

    # Check if hydrazide formation is in the late stage (close to the minimum depth)
    has_late_stage_hydrazide = False
    if hydrazide_formation_reaction:
        hydrazide_depth = hydrazide_formation_reaction[1]
        # Consider it late-stage if it's within 2 steps of the latest stage reaction
        has_late_stage_hydrazide = hydrazide_depth <= late_stage_depth + 2
        if has_late_stage_hydrazide:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "hydrazide_formation", "position": "late_stage", "description": "Hydrazide formation must occur within 2 steps of the latest reaction in the synthesis path (i.e., closest to the final product)."}})

    # For the strategy to be present, we need multiple sequential SNAr reactions and late-stage hydrazide formation
    strategy_present = has_multiple_snar and are_sequential and has_late_stage_hydrazide

    if strategy_present:
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["sequential_SNAr_reactions", "late_stage_hydrazide_formation"], "description": "The overall strategy is valid only if both a sequence of SNAr reactions and a late-stage hydrazide formation are present in the route."}})

    print(f"Sequential SNAr strategy detected: {strategy_present}")
    print(f"Number of SNAr reactions: {len(snar_reactions)}")
    print(f"SNAr reactions are sequential: {are_sequential}")
    print(f"Late-stage hydrazide formation: {has_late_stage_hydrazide}")

    return strategy_present, findings_json
