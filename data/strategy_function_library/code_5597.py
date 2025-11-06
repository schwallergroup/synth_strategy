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


HETEROCYCLES_OF_INTEREST = [
    "furan", "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole",
    "thiazole", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "indole", "quinoline", "isoquinoline", "benzimidazole", "benzoxazole",
    "benzothiazole", "triazole", "tetrazole", "oxadiazole", "thiadiazole",
    "pyrimidine", "pyrazine", "pyridazine", "aziridine", "azetidine",
    "pyrrolidine", "oxirane", "oxetane", "tetrahydrofuran", "tetrahydropyran",
    "dioxane", "thiazolidine", "oxazolidine", "isoxazole", "isothiazole"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving at least two N-alkylation steps, a late-stage fragment coupling reaction, and the formation of a specific heterocycle from a predefined list (e.g., pyrrole, pyridine, imidazole).
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

    n_alkylation_count = 0
    heterocycle_formation = False
    late_stage_fragment_coupling_with_linker = False

    def is_n_alkylation(reaction_smiles, reactants_smiles, product_smiles):
        """Check if the reaction is a known N-alkylation reaction."""
        n_alkylation_reactions = [
            "N-alkylation of primary amines with alkyl halides",
            "N-alkylation of secondary amines with alkyl halides",
            "Alkylation of amines",
            "Reductive amination with aldehyde",
            "Reductive amination with ketone",
            "N-methylation",
            "Eschweiler-Clarke Primary Amine Methylation",
            "Eschweiler-Clarke Secondary Amine Methylation",
            "Reductive methylation of primary amine with formaldehyde"
        ]
        for reaction_name in n_alkylation_reactions:
            if checker.check_reaction(reaction_name, reaction_smiles):
                if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                return True
        return False

    def is_heterocycle_formation(reactants_smiles, product_smiles, rsmi):
        """Check if the reaction forms a heterocycle."""
        for ring in HETEROCYCLES_OF_INTEREST:
            if checker.check_ring(ring, product_smiles):
                ring_in_reactants = False
                for r in reactants_smiles:
                    if checker.check_ring(ring, r):
                        ring_in_reactants = True
                        break
                if not ring_in_reactants:
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    return True

        heterocycle_reactions = [
            "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
            "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "Huisgen 1,3 dipolar cycloaddition",
            "Huisgen alkene-azide 1,3 dipolar cycloaddition", "Pyrazole formation",
            "Intramolecular amination (heterocycle formation)", "{benzimidazole_derivatives_carboxylic-acid/ester}",
            "{benzimidazole_derivatives_aldehyde}", "{benzothiazole}", "{benzoxazole_arom-aldehyde}",
            "{benzoxazole_carboxylic-acid}", "{thiazole}", "{tetrazole_terminal}",
            "{tetrazole_connect_regioisomere_1}", "{tetrazole_connect_regioisomere_2}",
            "{1,2,4-triazole_acetohydrazide}", "{1,2,4-triazole_carboxylic-acid/ester}",
            "{pyrazole}", "{Paal-Knorr pyrrole}", "{oxadiazole}", "{imidazole}",
        ]
        for reaction in heterocycle_reactions:
            if checker.check_reaction(reaction, rsmi):
                if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reaction)
                return True

        if checker.check_fg("Oxadiazole", product_smiles) and not any(
            checker.check_fg("Oxadiazole", r) for r in reactants_smiles
        ):
            if "Oxadiazole" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Oxadiazole")
            return True

        return False

    def has_linker_chain(reactants_smiles, product_smiles, rsmi):
        """Check if the reaction is a known coupling reaction."""
        coupling_reactions = [
            "Williamson Ether Synthesis", "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
            "Suzuki coupling with boronic acids", "Heck terminal vinyl",
            "Sonogashira acetylene_aryl halide", "{Williamson ether}", "{Suzuki}",
            "{Buchwald-Hartwig}", "{N-arylation_heterocycles}",
        ]
        for reaction in coupling_reactions:
            if checker.check_reaction(reaction, rsmi):
                if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reaction)
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_count, heterocycle_formation, late_stage_fragment_coupling_with_linker, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if is_n_alkylation(rsmi, reactants_smiles, product_smiles):
                    n_alkylation_count += 1

                if is_heterocycle_formation(reactants_smiles, product_smiles, rsmi):
                    heterocycle_formation = True

                if depth <= 2 and len(reactants_smiles) > 1:
                    if has_linker_chain(reactants_smiles, product_smiles, rsmi):
                        late_stage_fragment_coupling_with_linker = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    strategy_detected = (
        n_alkylation_count >= 2
        and heterocycle_formation
        and late_stage_fragment_coupling_with_linker
    )

    # Record structural constraints if detected
    if n_alkylation_count >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "N-alkylation",
                "operator": ">=",
                "value": 2
            }
        })
    
    if late_stage_fragment_coupling_with_linker:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "fragment_coupling",
                "position": "late_stage"
            }
        })

    if n_alkylation_count >= 2 and heterocycle_formation and late_stage_fragment_coupling_with_linker:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "N-alkylation",
                    "heterocycle_formation",
                    "late_stage_fragment_coupling"
                ]
            }
        })

    return strategy_detected, findings_json
