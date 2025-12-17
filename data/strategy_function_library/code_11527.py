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


NITROGEN_HETEROCYCLES = [
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "aziridine", "azetidine", "azepane", "diazepane", "indole", "quinoline",
    "isoquinoline", "purine", "carbazole", "acridine", "benzoxazole",
    "benzothiazole", "benzimidazole", "indazole", "benzotriazole",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzimidazole_derivatives_aldehyde}", "{benzothiazole}",
    "{benzoxazole_arom-aldehyde}", "{benzoxazole_carboxylic-acid}",
    "{thiazole}", "{tetrazole_terminal}", "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}", "{1,2,4-triazole_acetohydrazide}",
    "{1,2,4-triazole_carboxylic-acid/ester}", "{pyrazole}",
    "{Paal-Knorr pyrrole}", "{triaryl-imidazole}", "{Fischer indole}",
    "{indole}", "{oxadiazole}", "Benzothiazole formation from aldehyde",
    "Benzothiazole formation from acyl halide",
    "Benzothiazole formation from ester/carboxylic acid",
    "Benzoxazole formation from aldehyde",
    "Benzoxazole formation from acyl halide",
    "Benzoxazole formation from ester/carboxylic acid",
    "Benzoxazole formation (intramolecular)",
    "Benzimidazole formation from aldehyde",
    "Benzimidazole formation from acyl halide",
    "Benzimidazole formation from ester/carboxylic acid",
    "Paal-Knorr pyrrole synthesis", "Pyrazole formation",
    "Intramolecular amination (heterocycle formation)",
]

NITROGEN_MANIPULATION_REACTIONS = [
    "Reduction of nitro groups to amines", "Reductive amination with aldehyde",
    "Reductive amination with ketone", "Reductive amination with alcohol",
    "Reduction of nitrile to amine", "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines", "Reduction of tertiary amides to amines",
    "Azide to amine reduction (Staudinger)", "Primary amine to azide",
    "Amine to azide", "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Acylation of primary amines", "Acylation of secondary amines",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation", "N-methylation",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with primary amine to imide",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide", "Ester with ammonia to amide",
    "Ester with primary amine to amide", "Ester with secondary amine to amide",
    "Boc amine deprotection", "Boc amine protection", "Alkylation of amines",
    "Reductive amination", "{reductive amination}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a combined strategy involving both late-stage formation of a nitrogen-containing heterocycle and a nitrogen functional group interconversion. Late-stage is defined as the final third of the synthesis (e.g., steps 1-3 in a 9-step synthesis). Heterocycle formation is identified by checking for specific named reactions or the de novo appearance of a ring from the `NITROGEN_HETEROCYCLES` list. Nitrogen functional group manipulation is identified by checking for specific named reactions from the `NITROGEN_MANIPULATION_REACTIONS` list.
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

    # Track if we found each strategy
    found_late_heterocyclization = False
    found_nitrogen_manipulation = False

    # Track depth during traversal
    max_depth = 0

    # First pass to determine max depth only
    def count_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            # Depth only increases when going from chemical to reaction
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth += 1
            count_depth(child, new_depth)

    # Helper function to check if a reaction forms a nitrogen heterocycle
    def forms_nitrogen_heterocycle(reaction_node):
        nonlocal findings_json
        try:
            rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check specific heterocycle formation reactions
            for rxn in HETEROCYCLE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    print(f"Found heterocycle formation reaction: {rxn}")
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    return True

            # Check if product has a nitrogen heterocycle that reactants don't have
            for ring in NITROGEN_HETEROCYCLES:
                if checker.check_ring(ring, product):
                    # Check if any reactant already has this specific ring
                    reactants_have_ring = False
                    for reactant in reactants:
                        if checker.check_ring(ring, reactant):
                            reactants_have_ring = True
                            break
                    if not reactants_have_ring:
                        print(f"Found heterocycle formation: {ring}")
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        return True

            return False
        except Exception as e:
            print(f"Error checking heterocycle formation: {e}")
            return False

    # Helper function to check if a reaction involves nitrogen functional group manipulation
    def is_nitrogen_fg_manipulation(reaction_node):
        nonlocal findings_json
        try:
            rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]

            # Check specific nitrogen manipulation reactions
            for rxn in NITROGEN_MANIPULATION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    print(f"Found nitrogen manipulation reaction: {rxn}")
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    return True

            return False
        except Exception as e:
            print(f"Error checking nitrogen FG manipulation: {e}")
            return False

    # DFS traversal to analyze the synthesis route
    def dfs_traverse(node, depth=0):
        nonlocal found_late_heterocyclization, found_nitrogen_manipulation, findings_json

        # Check for nitrogen manipulation in reaction nodes
        if node["type"] == "reaction":
            # Check if this reaction forms a nitrogen heterocycle
            if forms_nitrogen_heterocycle(node):
                # Late stage is typically the first 1/3 of steps in retrosynthesis
                if depth <= max_depth / 3:
                    found_late_heterocyclization = True
                    print(f"Found late-stage heterocyclization at depth {depth}/{max_depth}")
                    # Add structural constraint for late_stage_nitrogen_heterocycle_formation
                    constraint_obj = {
                        "type": "positional",
                        "details": {
                            "target": "nitrogen_heterocycle_formation",
                            "position": "late_stage",
                            "condition": "depth <= max_depth / 3"
                        }
                    }
                    if constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_obj)

            # Check for nitrogen functional group manipulation
            if is_nitrogen_fg_manipulation(node):
                found_nitrogen_manipulation = True
                print(f"Found nitrogen functional group manipulation at depth {depth}/{max_depth}")

        # Recursively process children
        for child in node.get("children", []):
            # Depth only increases when going from chemical to reaction
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth += 1
            dfs_traverse(child, new_depth)

    # First determine the maximum depth
    count_depth(route)
    print(f"Maximum synthesis depth: {max_depth}")

    # Then analyze the route with known max depth
    dfs_traverse(route)

    # Combined strategy is present if both individual strategies are detected
    combined_strategy = found_late_heterocyclization and found_nitrogen_manipulation

    if combined_strategy:
        print("Detected combined strategy: late-stage heterocyclization with nitrogen manipulation")
        # Add structural constraint for co-occurrence
        constraint_obj = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "late_stage_nitrogen_heterocycle_formation",
                    "nitrogen_functional_group_interconversion"
                ]
            }
        }
        if constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint_obj)
    else:
        if found_late_heterocyclization:
            print("Found late-stage heterocyclization but no nitrogen manipulation")
        elif found_nitrogen_manipulation:
            print("Found nitrogen manipulation but no late-stage heterocyclization")
        else:
            print("Neither strategy was detected")

    return combined_strategy, findings_json
