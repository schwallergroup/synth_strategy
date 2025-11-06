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


LATE_STAGE_AMIDE_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with secondary amine to amide",
    "Ester with ammonia to amide",
    "Ester with secondary amine to amide",
    "Nitrile to amide",
    "Carboxylic acid to amide conversion",
]

AZIDE_REDUCTION_REACTIONS = [
    "Azide to amine reduction (Staudinger)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where an amine, formed from an azide intermediate, is used in a late-stage (final 3 steps) amide formation. This is confirmed by identifying an azide reduction reaction (specifically 'Azide to amine reduction (Staudinger)') followed by an amide formation reaction using the resulting amine. The function checks for a specific list of amide formation reactions: Acylation of Nitrogen Nucleophiles by Carboxylic Acids, Carboxylic acid with primary amine to amide, Ester with primary amine to amide, Acyl chloride with primary amine to amide (Schotten-Baumann), Schotten-Baumann_amide, Acylation of primary amines, Acylation of secondary amines, Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N, Acyl chloride with ammonia to amide, Acyl chloride with secondary amine to amide, Ester with ammonia to amide, Ester with secondary amine to amide, Nitrile to amide, Carboxylic acid to amide conversion. The connection between the amine intermediate and its use in the amide coupling is verified by SMILES or substructure matching.
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

    has_late_stage_amide = False
    has_azide_intermediate = False
    has_azide_to_amine_conversion = False

    azide_molecules = set()
    amine_from_azide = set()
    amines_used_in_amide = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_amide, has_azide_intermediate, has_azide_to_amine_conversion, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if depth <= 2:
                    for reaction_type in LATE_STAGE_AMIDE_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            has_late_stage_amide = True
                            if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                            if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            for reactant in reactants:
                                if checker.check_fg("Primary amine", reactant):
                                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                                    amines_used_in_amide.add(reactant)
                                elif checker.check_fg("Secondary amine", reactant):
                                    if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                                    amines_used_in_amide.add(reactant)

                    if not has_late_stage_amide:
                        if (checker.check_fg("Primary amide", product) or checker.check_fg("Secondary amide", product) or checker.check_fg("Tertiary amide", product)):
                            if checker.check_fg("Primary amide", product) and "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                            if checker.check_fg("Secondary amide", product) and "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                            if checker.check_fg("Tertiary amide", product) and "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")
                            for reactant in reactants:
                                if checker.check_fg("Primary amine", reactant):
                                    has_late_stage_amide = True
                                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                                    amines_used_in_amide.add(reactant)
                                elif checker.check_fg("Secondary amine", reactant):
                                    has_late_stage_amide = True
                                    if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                                    amines_used_in_amide.add(reactant)
                            if has_late_stage_amide and "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                            if has_late_stage_amide and {"type": "positional", "details": {"target": "amide_formation", "position": "last_3_stages", "description": "The amide formation reaction must occur within the last three steps of the synthesis (depth <= 2)."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "last_3_stages", "description": "The amide formation reaction must occur within the last three steps of the synthesis (depth <= 2)."}})

                for smiles in [product] + reactants:
                    if checker.check_fg("Azide", smiles):
                        has_azide_intermediate = True
                        if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Azide")
                        azide_molecules.add(smiles)

                for reaction_type in AZIDE_REDUCTION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        if any(checker.check_fg("Azide", r) for r in reactants) and (checker.check_fg("Primary amine", product) or checker.check_fg("Secondary amine", product)):
                            has_azide_to_amine_conversion = True
                            if "azide_reduction" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("azide_reduction")
                            if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            if checker.check_fg("Primary amine", product) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            if checker.check_fg("Secondary amine", product) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            amine_from_azide.add(product)

                if not has_azide_to_amine_conversion:
                    if (checker.check_fg("Primary amine", product) or checker.check_fg("Secondary amine", product)) and any(checker.check_fg("Azide", r) for r in reactants):
                        has_azide_to_amine_conversion = True
                        if "azide_reduction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("azide_reduction")
                        if checker.check_fg("Primary amine", product) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", product) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        amine_from_azide.add(product)

        elif node["type"] == "mol":
            smiles = node["smiles"]
            if checker.check_fg("Azide", smiles):
                has_azide_intermediate = True
                if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")
                azide_molecules.add(smiles)

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    amine_azide_connection = False
    if amine_from_azide.intersection(amines_used_in_amide):
        amine_azide_connection = True
        if {"type": "sequence", "details": {"first_event": "azide_reduction", "second_event": "amide_formation", "linkage": "product_is_reactant", "description": "The amine product of the azide reduction must be a reactant in the subsequent amide formation."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"first_event": "azide_reduction", "second_event": "amide_formation", "linkage": "product_is_reactant", "description": "The amine product of the azide reduction must be a reactant in the subsequent amide formation."}})

    if not amine_azide_connection and amine_from_azide and amines_used_in_amide:
        try:
            for amine_azide in amine_from_azide:
                for amine_amide in amines_used_in_amide:
                    mol1 = Chem.MolFromSmiles(amine_azide)
                    mol2 = Chem.MolFromSmiles(amine_amide)
                    if mol1 and mol2:
                        mcs = rdFMCS.FindMCS([mol1, mol2], bondCompare=rdFMCS.BondCompare.CompareOrder, atomCompare=rdFMCS.AtomCompare.CompareElements, ringMatchesRingOnly=True, completeRingsOnly=True)
                        if mcs.numAtoms > 0:
                            smaller_atom_count = min(mol1.GetNumAtoms(), mol2.GetNumAtoms())
                            if mcs.numAtoms >= 0.7 * smaller_atom_count:
                                amine_azide_connection = True
                                if {"type": "sequence", "details": {"first_event": "azide_reduction", "second_event": "amide_formation", "linkage": "product_is_reactant", "description": "The amine product of the azide reduction must be a reactant in the subsequent amide formation."}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"first_event": "azide_reduction", "second_event": "amide_formation", "linkage": "product_is_reactant", "description": "The amine product of the azide reduction must be a reactant in the subsequent amide formation."}})
                                break
                if amine_azide_connection:
                    break
        except Exception:
            pass

    result = (
        has_late_stage_amide
        and has_azide_intermediate
        and has_azide_to_amine_conversion
        and amine_azide_connection
    )

    if result:
        if {"type": "co-occurrence", "details": {"targets": ["amide_formation", "azide_reduction", "Azide"], "description": "The route must contain an amide formation, an azide reduction, and an azide functional group."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["amide_formation", "azide_reduction", "Azide"], "description": "The route must contain an amide formation, an azide reduction, and an azide functional group."}})

    return result, findings_json
