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


string:HETEROCYCLE_TYPES_FOR_COUPLING = [
    "furan", "pyrrole", "thiophene", "pyrazole", "imidazole", "oxazole",
    "thiazole", "triazole", "tetrazole", "pyridine", "pyrimidine", "pyrazine",
    "pyridazine", "indole", "benzimidazole", "benzothiazole", "benzoxazole",
    "quinoline", "isoquinoline", "purine", "morpholine", "piperidine",
    "piperazine", "pyrrolidine", "oxazolidine", "thiazolidine", "isoxazole",
    "isothiazole",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Ester with primary amine to amide",
    "Acyl chloride with secondary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Schotten-Baumann to ester",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage coupling of two distinct heterocyclic fragments via amide bond formation. The specific reaction types and heterocycles of interest are defined in the `AMIDE_FORMATION_REACTIONS` and `HETEROCYCLE_TYPES_FOR_COUPLING` lists, respectively. The check is limited to the final two synthetic steps.
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

    found_heterocycle_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_coupling, findings_json

        # Only check reaction nodes in the final two steps (depth=1 and depth=2)
        if node["type"] == "reaction" and depth <= 2:
            # Record positional constraint if met
            if depth <= 2:
                if {"type": "positional", "details": {"target": "amide_coupling_of_heterocycles", "position": "last_two_stages"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_coupling_of_heterocycles", "position": "last_two_stages"}})

            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                is_amide_formation = False
                # Check if this is an amide formation reaction by name
                for reaction_type in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amide_formation = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                # Fallback: Check if reaction forms an amide bond by FG analysis
                if not is_amide_formation:
                    has_amide_product = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )
                    if has_amide_product:
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Primary amide", product_smiles):
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Secondary amide", product_smiles):
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"] and checker.check_fg("Tertiary amide", product_smiles):
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    has_amide_reactants = False
                    for r in reactants_smiles:
                        if checker.check_fg("Primary amide", r):
                            has_amide_reactants = True
                            if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        if checker.check_fg("Secondary amide", r):
                            has_amide_reactants = True
                            if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        if checker.check_fg("Tertiary amide", r):
                            has_amide_reactants = True
                            if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    has_acid_or_ester = False
                    for r in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", r):
                            has_acid_or_ester = True
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if checker.check_fg("Ester", r):
                            has_acid_or_ester = True
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        if checker.check_fg("Acyl halide", r):
                            has_acid_or_ester = True
                            if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

                    has_amine = False
                    for r in reactants_smiles:
                        if checker.check_fg("Primary amine", r):
                            has_amine = True
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", r):
                            has_amine = True
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if checker.check_fg("Aniline", r):
                            has_amine = True
                            if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                    if (has_amide_product and not has_amide_reactants and has_acid_or_ester and has_amine):
                        is_amide_formation = True
                        # Record structural constraints for amide formation by FG analysis
                        if {"type": "co-occurrence", "details": {"targets": ["amide_in_product", "acid_or_ester_in_reactant", "amine_in_reactant"]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["amide_in_product", "acid_or_ester_in_reactant", "amine_in_reactant"]}})
                        if {"type": "negation", "details": {"target": "amide_in_reactant"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "amide_in_reactant"}})
                        if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("amide_formation")

                if not is_amide_formation:
                    return

                # Check for heterocycles in reactants
                heterocyclic_reactants = []
                for reactant in reactants_smiles:
                    for heterocycle in HETEROCYCLE_TYPES_FOR_COUPLING:
                        if checker.check_ring(heterocycle, reactant):
                            heterocyclic_reactants.append(reactant)
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            break

                # A true coupling requires at least two different reactant molecules with heterocycles.
                # Using a set of reactant SMILES ensures we count distinct fragments.
                if len(set(heterocyclic_reactants)) >= 2:
                    found_heterocycle_coupling = True
                    # Record structural constraint for distinct heterocyclic reactants
                    if {"type": "count", "details": {"target": "distinct_heterocyclic_reactants", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "distinct_heterocyclic_reactants", "operator": ">=", "value": 2}})

            except Exception:
                # Silently ignore errors in processing a single node
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return found_heterocycle_coupling, findings_json
