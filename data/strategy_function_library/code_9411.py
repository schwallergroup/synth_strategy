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


# Refactored module-level constants
HETEROCYCLES_OF_INTEREST = [
    "benzimidazole", "imidazole", "oxazole", "thiazole", "pyrazole",
    "isoxazole", "isothiazole", "triazole", "tetrazole", "pyridine",
    "pyrimidine", "pyrazine", "pyridazine", "indole", "benzoxazole",
    "benzothiazole", "furan", "thiophene", "pyrrole", "oxadiazole",
    "thiadiazole", "morpholine", "piperidine", "piperazine",
]

AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann to ester",
    "{Schotten-Baumann_amide}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a common synthetic pattern involving the early-stage formation of a specific heterocycle
    followed by a late-stage amide coupling. The specific heterocycles and amide coupling reaction
    names are defined in the HETEROCYCLES_OF_INTEREST and AMIDE_COUPLING_REACTIONS lists, respectively.
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

    heterocycle_formed_early = False
    amide_coupling_late = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed_early, amide_coupling_late, findings_json

        if node["type"] == "reaction":
            # Extract product and reactants
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocycle formation
            product_has_heterocycle = False
            heterocycle_found = None

            for heterocycle in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(heterocycle, product_smiles):
                    product_has_heterocycle = True
                    heterocycle_found = heterocycle
                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                    break

            if product_has_heterocycle:
                # Check if any reactant already has the heterocycle
                reactants_have_heterocycle = False
                for reactant_smiles in reactants_smiles:
                    if checker.check_ring(heterocycle_found, reactant_smiles):
                        reactants_have_heterocycle = True
                        break

                # Early stage is defined as depth > 1
                if not reactants_have_heterocycle and depth > 1:
                    heterocycle_formed_early = True
                    # Record positional constraint for ring formation
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "ring_formation",
                            "position": "early_stage"
                        }
                    })

            # Check for amide coupling reactions
            for reaction_type in AMIDE_COUPLING_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    # Late stage is defined as depth <= 1
                    if depth <= 1:
                        amide_coupling_late = True
                        # Record positional constraint for amide coupling
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "amide_coupling",
                                "position": "late_stage"
                            }
                        })
                    break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    combined_strategy = heterocycle_formed_early and amide_coupling_late
    if combined_strategy:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "early_stage_heterocycle_formation",
                    "late_stage_amide_coupling"
                ]
            }
        })

    return combined_strategy, findings_json