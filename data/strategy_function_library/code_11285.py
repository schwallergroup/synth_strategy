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

HETEROCYCLES_OF_INTEREST = [
    "thiophene", "furan", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "benzoxazole", "benzothiazole", "benzimidazole", "indole",
    "quinoline", "isoquinoline", "piperidine", "piperazine", "morpholine",
    "thiomorpholine", "oxadiazole", "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy that involves both the formation of a key heterocycle and a late-stage amide coupling.
    This pattern is common in convergent syntheses where a heterocyclic core is built and then coupled to another fragment.
    This check identifies the formation of one of the following heterocycles: thiophene, furan, pyrrole, pyridine, pyrazole, imidazole, oxazole, thiazole, pyrimidine, pyrazine, pyridazine, triazole, tetrazole, benzoxazole, benzothiazole, benzimidazole, indole, quinoline, isoquinoline, piperidine, piperazine, morpholine, thiomorpholine, oxadiazole, or thiadiazole.
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

    amide_coupling_at_late_stage = False
    heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_at_late_stage, heterocycle_formation, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for amide coupling at late stage (depth 0-3)
                if depth <= 3:
                    amide_coupling_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with primary amine to amide",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with secondary amine to amide",
                        "Acyl chloride with ammonia to amide",
                        "Ester with ammonia to amide",
                        "Schotten-Baumann_amide",
                        "Acylation of primary amines",
                        "Acylation of secondary amines"
                    ]
                    is_amide_coupling = False
                    for r_name in amide_coupling_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            is_amide_coupling = True
                            if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            break

                    if is_amide_coupling:
                        amide_coupling_at_late_stage = True

                # Check for heterocycle formation at early to mid stage
                if depth >= 2:
                    product_has_heterocycle = False
                    for ring in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring, product_smiles):
                            product_has_heterocycle = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                    reactants_have_heterocycle = False
                    for r in reactants_smiles:
                        for ring in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring, r):
                                reactants_have_heterocycle = True
                                # Do not add to findings_json here, as it's about formation, not presence in reactants
                                break
                        if reactants_have_heterocycle: break

                    heterocycle_formation_reactions = [
                        "Formation of NOS Heterocycles",
                        "Paal-Knorr pyrrole synthesis",
                        "benzimidazole_derivatives_aldehyde",
                        "benzimidazole_derivatives_carboxylic-acid/ester",
                        "benzothiazole",
                        "benzoxazole_arom-aldehyde",
                        "benzoxazole_carboxylic-acid",
                        "thiazole",
                        "tetrazole_terminal",
                        "tetrazole_connect_regioisomere_1",
                        "tetrazole_connect_regioisomere_2",
                        "1,2,4-triazole_acetohydrazide",
                        "1,2,4-triazole_carboxylic-acid/ester",
                        "pyrazole",
                        "Fischer indole",
                        "oxadiazole",
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        "Huisgen 1,3 dipolar cycloaddition",
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                        "Pyrazole formation"
                    ]
                    is_heterocycle_formation_reaction = False
                    for r_name in heterocycle_formation_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            is_heterocycle_formation_reaction = True
                            if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            break

                    if (product_has_heterocycle and not reactants_have_heterocycle) or is_heterocycle_formation_reaction:
                        heterocycle_formation = True

            except Exception as e:
                # In a production setting, you might want to log this error.
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = amide_coupling_at_late_stage and heterocycle_formation

    if result:
        # Add the structural constraint if both conditions are met
        # This corresponds to the overall strategy being detected
        structural_constraint_obj = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "amide_coupling",
                    "heterocycle_formation"
                ],
                "constraints": [
                    {
                        "type": "positional",
                        "details": {
                            "target": "amide_coupling",
                            "position": "late_stage",
                            "condition": "depth <= 3"
                        }
                    },
                    {
                        "type": "positional",
                        "details": {
                            "target": "heterocycle_formation",
                            "position": "early_to_mid_stage",
                            "condition": "depth >= 2"
                        }
                    }
                ]
            }
        }
        if structural_constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(structural_constraint_obj)

    return result, findings_json
