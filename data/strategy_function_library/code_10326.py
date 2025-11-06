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


# Refactored lists as module-level constants
HETEROCYCLE_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "benzimidazole_derivatives_carboxylic-acid/ester", "benzimidazole_derivatives_aldehyde",
    "benzothiazole", "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid",
    "thiazole", "Niementowski_quinazoline", "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1", "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide", "1,2,4-triazole_carboxylic-acid/ester",
    "3-nitrile-pyridine", "pyrazole", "phthalazinone", "Fischer indole",
    "Friedlaender chinoline", "benzofuran", "benzothiophene", "indole",
    "oxadiazole", "imidazole", "Pictet-Spengler",
]

REDUCTIVE_AMINATION_REACTIONS = [
    "reductive amination with aldehyde", "reductive amination with ketone",
    "reductive amination with alcohol", "reductive amination",
]

HETEROCYCLES_FOR_FORMATION_CHECK = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "trioxane", "dioxepane", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "pyrrolidine", "piperidine", "piperazine", "morpholine",
    "thiomorpholine", "aziridine", "azetidine", "azepane", "diazepane",
    "indole", "quinoline", "isoquinoline", "purine", "carbazole", "acridine",
    "thiophene", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dithiane", "dithiolane", "benzothiophene", "oxathiolane", "dioxathiolane",
    "thiazolidine", "oxazolidine", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "benzene", "naphthalene", "anthracene", "benzoxazole",
    "benzothiazole", "benzimidazole", "pteridin", "phenothiazine",
    "phenoxazine", "dibenzofuran", "dibenzothiophene", "xanthene",
    "thioxanthene", "pyrroline", "pyrrolidone", "imidazolidine", "porphyrin",
    "indazole", "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a reductive amination occurs at a later stage
    than a heterocycle formation step. Reductive amination is identified from the
    list `REDUCTIVE_AMINATION_REACTIONS`. Heterocycle formation is identified either
    from the list `HETEROCYCLE_FORMATION_REACTIONS` or by detecting the de novo
    formation of a ring specified in `HETEROCYCLES_FOR_FORMATION_CHECK`.
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

    all_reactions = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal all_reactions, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reaction_info = {"depth": depth, "type": "unknown"}

                # Check for reductive amination by name
                for rxn in REDUCTIVE_AMINATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        reaction_info["type"] = "reductive_amination"
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break
                
                # If not reductive amination, check for heterocycle formation
                if reaction_info["type"] == "unknown":
                    # Check by reaction name
                    for rxn in HETEROCYCLE_FORMATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            reaction_info["type"] = "heterocycle_formation"
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break
                
                if reaction_info["type"] == "unknown":
                    # Check by de novo ring formation
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]
                    for ring in HETEROCYCLES_FOR_FORMATION_CHECK:
                        is_formed = checker.check_ring(ring, product_smiles) and \
                                    not any(checker.check_ring(ring, r) for r in reactants_smiles)
                        if is_formed:
                            reaction_info["type"] = "heterocycle_formation"
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            break
                
                if reaction_info["type"] != "unknown":
                    all_reactions.append(reaction_info)

            except Exception:
                # Silently ignore errors in processing a single reaction
                pass

        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    reductive_aminations = [r for r in all_reactions if r["type"] == "reductive_amination"]
    heterocycle_formations = [r for r in all_reactions if r["type"] == "heterocycle_formation"]

    if not reductive_aminations or not heterocycle_formations:
        return False, findings_json

    # Strategy is true if ANY reductive amination is later than ANY heterocycle formation
    # Later step in forward synthesis means smaller depth in retrosynthesis
    for ra in reductive_aminations:
        for hf in heterocycle_formations:
            if ra["depth"] < hf["depth"]:
                result = True
                # Add structural constraints if the condition is met
                if {"type": "co-occurrence", "details": {"targets": ["reductive_amination", "heterocycle_formation"]}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["reductive_amination", "heterocycle_formation"]}})
                if {"type": "sequence", "details": {"before": "heterocycle_formation", "after": "reductive_amination"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "heterocycle_formation", "after": "reductive_amination"}})
                return result, findings_json

    return result, findings_json