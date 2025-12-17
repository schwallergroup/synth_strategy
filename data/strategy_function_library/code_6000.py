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
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "imidazole", "oxazole", "thiazole", "furan", "thiophene", "pyrrole", "indole",
    "benzimidazole", "benzoxazole", "benzothiazole", "quinoline", "isoquinoline",
    "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran", "oxirane", "oxetane",
    "oxolane", "oxane", "dioxolane", "dioxolene", "trioxane", "dioxepane", "pyrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "aziridine", "azetidine", "azepane", "diazepane", "purine", "carbazole",
    "acridine", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dithiane", "dithiolane", "benzothiophene", "oxathiolane", "dioxathiolane",
    "thiazolidine", "oxazolidine", "isoxazole", "isothiazole", "oxadiazole", "thiadiazole"
]

BIARYL_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Stille reaction_aryl", "Stille reaction_vinyl", "Negishi coupling",
    "Kumada cross-coupling", "Hiyama-Denmark Coupling", "Ullmann condensation",
    "Suzuki coupling with boronic acids OTf", "Suzuki coupling with sulfonic esters",
    "Suzuki coupling with boronic esters OTf", "Aryllithium cross-coupling"
]

FUNCTIONALIZATION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides", "N-alkylation of secondary amines with alkyl halides",
    "Aromatic bromination", "Aromatic chlorination", "Aromatic iodination", "Aromatic fluorination",
    "N-arylation", "Buchwald-Hartwig", "Acylation of Nitrogen Nucleophiles",
    "Friedel-Crafts acylation", "Friedel-Crafts alkylation", "Aromatic nitration",
    "Methylation", "Alkylation of amines", "Acylation of primary amines", "Acylation of secondary amines",
    "Esterification of Carboxylic Acids", "Williamson Ether Synthesis", "Reductive amination with aldehyde",
    "Reductive amination with ketone", "Oxidation of aldehydes to carboxylic acids",
    "Oxidation of alcohol to carboxylic acid", "Oxidation of alcohol to aldehyde",
    "Oxidation of ketone to carboxylic acid", "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol", "Reduction of carboxylic acid to primary alcohol",
    "Reduction of nitrile to amine", "Reduction of nitro groups to amines", "Minisci (para)",
    "Minisci (ortho)", "Directed ortho metalation of arenes", "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt", "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2", "Aromatic hydroxylation", "Methylation with MeI_primary",
    "Methylation with MeI_secondary", "Methylation with MeI_tertiary", "Methylation with MeI_aryl",
    "Methylation with MeI_SH", "Methylation with DMS", "Methylation with DMC"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a 'heterocycle-first' strategy involving early-stage heterocycle formation, subsequent functionalization, and a late-stage biaryl coupling. The strategy is identified by checking for specific, enumerated lists of heterocycles (HETEROCYCLES_OF_INTEREST), functionalization reactions (FUNCTIONALIZATION_REACTIONS), and biaryl coupling reactions (BIARYL_COUPLING_REACTIONS).
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

    # Track key features
    heterocycle_formation_early = False
    heterocycle_present_in_starting = False
    late_biaryl_coupling = False
    sequential_functionalization = 0

    # Track reaction depths and heterocycles
    reaction_depths = {}
    max_depth = 0
    detected_heterocycles = set()

    # First pass: determine max depth and check for heterocycles in starting materials
    def first_pass(node, depth=0):
        nonlocal max_depth, heterocycle_present_in_starting, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            # Check if molecule contains heterocycle
            for heterocycle in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(heterocycle, node["smiles"]):
                    detected_heterocycles.add(heterocycle)
                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                    if node.get("in_stock", False):
                        heterocycle_present_in_starting = True

        # Process children
        for child in node.get("children", []):
            first_pass(child, depth + 1)

    # Main traversal
    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_early, late_biaryl_coupling, sequential_functionalization, findings_json

        if node["type"] == "reaction":
            # Store reaction at this depth
            reaction_depths[depth] = node

            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for heterocycle formation
                product_heterocycles = set()
                reactant_heterocycles = set()

                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_heterocycles.add(heterocycle)
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                    for reactant in reactants_smiles:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.add(heterocycle)
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                # Check for heterocycle formation (early stage)
                if depth >= max_depth / 2:
                    new_heterocycles = product_heterocycles - reactant_heterocycles
                    if new_heterocycles:
                        heterocycle_formation_early = True
                        detected_heterocycles.update(new_heterocycles)
                        # Add 'ring_formation' to named_reactions if a new heterocycle is formed
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # Check for biaryl coupling (late stage)
                if depth <= 2:  # Final, penultimate, or antepenultimate step
                    for coupling in BIARYL_COUPLING_REACTIONS:
                        if checker.check_reaction(coupling, rsmi):
                            late_biaryl_coupling = True
                            if coupling not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(coupling)
                            break

                # Check for functionalization steps
                if product_heterocycles and (
                    depth < max_depth / 2
                    or detected_heterocycles.intersection(product_heterocycles)
                ):
                    for reaction in FUNCTIONALIZATION_REACTIONS:
                        if checker.check_reaction(reaction, rsmi):
                            sequential_functionalization += 1
                            if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction)
                            break

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Run first pass to determine max depth and check starting materials
    first_pass(route)

    # Start main traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        (heterocycle_formation_early or heterocycle_present_in_starting)
        and late_biaryl_coupling
        and sequential_functionalization >= 1
    )

    # If we have a heterocycle and biaryl coupling but no functionalization,
    # check if there are any reactions between heterocycle formation and biaryl coupling
    if (
        (heterocycle_formation_early or heterocycle_present_in_starting)
        and late_biaryl_coupling
        and sequential_functionalization == 0
    ):
        if len(reaction_depths) > 2:
            sequential_functionalization = 1
            strategy_present = True

    # Populate structural constraints based on detected flags
    if heterocycle_formation_early or heterocycle_present_in_starting:
        # This corresponds to 'heterocycle_source'
        # If heterocycle_formation_early is true, it also implies 'ring_formation' and 'early_stage' positional constraint
        if heterocycle_formation_early:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "ring_formation",
                    "position": "early_stage",
                    "comment": "A possible source for the heterocycle is its de novo formation in the first half of the synthesis (depth >= max_depth / 2)."
                }
            })

    if late_biaryl_coupling:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target_group": "BIARYL_COUPLING_REACTIONS",
                "position": "late_stage",
                "comment": "A biaryl coupling reaction must occur within the last three steps of the synthesis (depth <= 2)."
            }
        })

    if sequential_functionalization >= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target_group": "FUNCTIONALIZATION_REACTIONS",
                "operator": ">=",
                "value": 1,
                "comment": "At least one reaction from the functionalization list must occur."
            }
        })

    if strategy_present:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "heterocycle_source",
                    "late_stage_biaryl_coupling",
                    "intermediate_functionalization"
                ],
                "comment": "The overall strategy requires a source for a heterocycle (either formed early or from starting materials), at least one functionalization reaction, and a late-stage biaryl coupling."
            }
        })

    return strategy_present, findings_json
