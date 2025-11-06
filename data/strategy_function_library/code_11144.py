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


HETEROCYCLES_FOR_RING_OPENING = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "pyrrole", "pyridine",
    "pyrazole", "imidazole", "oxazole", "thiazole", "pyrimidine", "pyrazine",
    "pyridazine", "triazole", "tetrazole", "pyrrolidine", "piperidine",
    "piperazine", "morpholine", "thiomorpholine", "aziridine", "azetidine",
    "thiophene", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dioxolane", "dioxolene", "trioxane", "dioxepane", "azepane", "diazepane",
    "indole", "quinoline", "isoquinoline", "benzothiophene", "oxathiolane",
    "dioxathiolane", "thiazolidine", "oxazolidine",
]

NUCLEOPHILIC_SUBSTITUTION_REACTIONS = [
    "Williamson Ether Synthesis", "S-alkylation of thiols",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Alcohol to azide", "thioether_nucl_sub", "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro", "nucl_sub_aromatic_para_nitro",
    "Mitsunobu aryl ether", "Mitsunobu esterification",
    "Mitsunobu aryl ether (intramolecular)", "Mitsunobu_imide",
    "Mitsunobu_phenole", "Mitsunobu_sulfonamide", "Mitsunobu_tetrazole_1",
    "Mitsunobu_tetrazole_2", "Mitsunobu_tetrazole_3", "Mitsunobu_tetrazole_4",
    "Ring opening of epoxide with amine", "Finkelstein reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a combined strategy where an early-stage ring opening of a heterocycle is followed by a late-stage nucleophilic substitution. Ring openings are identified when a heterocycle from the predefined `HETEROCYCLES_FOR_RING_OPENING` list is present in a reactant but not the product. Nucleophilic substitutions are identified by matching against the predefined `NUCLEOPHILIC_SUBSTITUTION_REACTIONS` list. An 'early-stage' reaction is defined as occurring at a depth greater than 2, and a 'late-stage' reaction at a depth of 3 or less.
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

    found_late_stage_substitution = False
    found_early_stage_ring_opening = False

    def is_nucleophilic_substitution(rsmi):
        nonlocal findings_json
        for rxn_type in NUCLEOPHILIC_SUBSTITUTION_REACTIONS:
            if checker.check_reaction(rxn_type, rsmi):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def is_ring_opening(reactants, product):
        nonlocal findings_json
        for ring in HETEROCYCLES_FOR_RING_OPENING:
            ring_found_in_reactant = False
            for reactant in reactants:
                if checker.check_ring(ring, reactant):
                    ring_found_in_reactant = True
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    break # Found in at least one reactant
            
            if ring_found_in_reactant:
                if not checker.check_ring(ring, product):
                    # This ring was in reactant but not product, indicating opening
                    # The ring itself is already added to findings_json["atomic_checks"]["ring_systems"]
                    return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_substitution, found_early_stage_ring_opening, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if is_nucleophilic_substitution(rsmi):
                    if depth <= 3:
                        found_late_stage_substitution = True
                        # Add positional constraint for late-stage nucleophilic substitution
                        constraint = {
                            "type": "positional",
                            "details": {
                                "target": "nucleophilic_substitution",
                                "position": "late_stage (depth <= 3)"
                            }
                        }
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)

                if is_ring_opening(reactants, product):
                    if depth > 2:
                        found_early_stage_ring_opening = True
                        # Add positional constraint for early-stage ring destruction
                        constraint = {
                            "type": "positional",
                            "details": {
                                "target": "ring_destruction",
                                "position": "early_stage (depth > 2)"
                            }
                        }
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)

            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = found_late_stage_substitution and found_early_stage_ring_opening
    
    if result:
        # Add co-occurrence constraint if both conditions are met
        constraint = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "late_stage_nucleophilic_substitution",
                    "early_stage_ring_destruction"
                ]
            }
        }
        if constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint)

    return result, findings_json
