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


HETEROCYCLE_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "tetrazole_terminal",
    "pyrazole",
    "imidazole",
    "oxadiazole",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "3-nitrile-pyridine",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen_disubst-alkyne",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Pictet-Spengler",
]

HETEROCYCLIC_RINGS_OF_INTEREST = [
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
    "thiadiazole", "benzoxazole", "benzothiazole", "benzimidazole", "pteridin",
    "phenothiazine", "phenoxazine", "dibenzofuran", "dibenzothiophene",
    "xanthene", "thioxanthene", "pyrroline", "pyrrolidone", "imidazolidine",
    "porphyrin", "indazole", "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage (final step) formation of a heterocycle. The strategy is identified if the reaction matches a known heterocycle-forming named reaction, or if a de novo formation of a heterocycle from an enumerated list is observed by comparing reactant and product structures.
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

    # Track if we found a heterocycle formation in the final step
    heterocycle_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found, findings_json

        if node["type"] == "reaction" and depth <= 1:
            # Extract reaction information
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Parse reactants and products
            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[-1]

            # Check if this is a heterocycle formation reaction
            for rxn_type in HETEROCYCLE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    heterocycle_formation_found = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    # Add the structural constraint if the condition is met at the final stage
                    if depth == 0: # Assuming depth 0 means final stage
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "heterocycle_formation",
                                "position": "last_stage"
                            }
                        })
                    return

            # If no specific reaction type was found, check for heterocycle formation by comparing rings
            try:
                # Check for rings in product that aren't in reactants
                for ring_name in HETEROCYCLIC_RINGS_OF_INTEREST:
                    if checker.check_ring(ring_name, product):
                        # Check if this ring wasn't in any reactant
                        ring_is_new = True
                        for reactant in reactants:
                            if checker.check_ring(ring_name, reactant):
                                ring_is_new = False
                                break

                        if ring_is_new:
                            heterocycle_formation_found = True
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            # Add the structural constraint if the condition is met at the final stage
                            if depth == 0: # Assuming depth 0 means final stage
                                findings_json["structural_constraints"].append({
                                    "type": "positional",
                                    "details": {
                                        "target": "heterocycle_formation",
                                        "position": "last_stage"
                                    }
                                })
                            return
            except Exception as e:
                pass

        # Process children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return heterocycle_formation_found, findings_json
