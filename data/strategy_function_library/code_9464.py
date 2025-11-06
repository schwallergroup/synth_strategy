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


# Module-level constants for heterocycle functionalization strategy
HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazine",
    "imidazole", "pyrazole", "oxazole", "thiazole", "triazole", "tetrazole",
    "indole", "benzimidazole", "purine", "quinoline", "isoquinoline",
    "furan", "thiophene", "pyrrole", "oxadiazole", "thiadiazole",
]

ACTIVATION_REACTIONS = [
    "Aromatic bromination", "Aromatic chlorination", "Aromatic iodination", "Aromatic fluorination",
    "Bromination", "Chlorination", "Iodination", "Fluorination",
    "Aromatic nitration with HNO3", "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt", "Aromatic nitration with alkyl NO2",
]

HALIDE_FGS_FOR_ACTIVATION = ["Primary halide", "Secondary halide", "Tertiary halide", "Aromatic halide"]

COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Sonogashira acetylene_aryl halide", "Sonogashira alkyne_aryl halide",
    "Heck terminal vinyl", "Negishi coupling", "Stille reaction_aryl",
    "Ullmann condensation", "Goldberg coupling", "Hiyama-Denmark Coupling",
    "Kumada cross-coupling", "Aryllithium cross-coupling",
    "Catellani reaction ortho", "Catellani reaction para",
    "Heck reaction with vinyl ester and amine", "Oxidative Heck reaction",
    "Oxidative Heck reaction with vinyl ester",
    "Reductive amination with aldehyde", "Reductive amination with ketone",
    "Reductive amination with alcohol",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Alkylation of amines", "Acylation of primary amines", "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
]

def main(route) -> Tuple[bool, Dict]:
    """Detects a two-step strategy for functionalizing a heterocyclic scaffold. The strategy consists of an initial activation step followed by a subsequent coupling reaction. The specific heterocycles, activation reactions, and coupling reactions checked are defined in the HETEROCYCLES_OF_INTEREST, ACTIVATION_REACTIONS, and COUPLING_REACTIONS lists, respectively. The function also handles cases where a pre-activated heterocycle (e.g., a halo-heterocycle) is used as a starting material."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Strategy detection flags
    has_heterocycle = False
    has_electrophilic_activation = False
    has_nucleophilic_coupling = False
    reaction_sequence_correct = False

    # Track reaction depths for sequence analysis
    activation_depth = -1
    coupling_depth = -1

    # Track pre-activated heterocycles
    pre_activated_heterocycles = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle, has_electrophilic_activation, has_nucleophilic_coupling
        nonlocal activation_depth, coupling_depth, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            for h in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(h, mol_smiles):
                    has_heterocycle = True
                    if h not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(h)

            # Check for pre-activated heterocycles (already containing halides)
            for h_fg in HALIDE_FGS_FOR_ACTIVATION:
                if checker.check_fg(h_fg, mol_smiles):
                    if h_fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(h_fg)
                    pre_activated_heterocycles.add(mol_smiles)
                    has_electrophilic_activation = True

        elif node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    is_heterocycle_reaction = any(
                        any(checker.check_ring(h, mol) for h in HETEROCYCLES_OF_INTEREST)
                        for mol in reactants + [product]
                    )

                    if is_heterocycle_reaction:
                        # Check for electrophilic activation reaction
                        if activation_depth == -1:
                            for r in ACTIVATION_REACTIONS:
                                if checker.check_reaction(r, rsmi):
                                    has_electrophilic_activation = True
                                    activation_depth = depth
                                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append(r)
                                    break
                            else:
                                # Fallback: check for halide formation
                                product_has_halide = any(checker.check_fg(h, product) for h in HALIDE_FGS_FOR_ACTIVATION)
                                reactants_have_halide = any(any(checker.check_fg(h, r) for h in HALIDE_FGS_FOR_ACTIVATION) for r in reactants)
                                if product_has_halide and not reactants_have_halide:
                                    has_electrophilic_activation = True
                                    activation_depth = depth
                                    if "halide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append("halide_formation")
                                    for h_fg in HALIDE_FGS_FOR_ACTIVATION:
                                        if checker.check_fg(h_fg, product) and h_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                            findings_json["atomic_checks"]["functional_groups"].append(h_fg)
                        
                        # Check for nucleophilic coupling reaction
                        if coupling_depth == -1:
                            for r in COUPLING_REACTIONS:
                                if checker.check_reaction(r, rsmi):
                                    has_nucleophilic_coupling = True
                                    coupling_depth = depth
                                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append(r)
                                    break
                            else:
                                # Fallback: check for N-functionalization
                                reactant_has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                                product_has_new_N_subst = checker.check_fg("Secondary amine", product) or checker.check_fg("Tertiary amine", product)
                                if reactant_has_primary_amine and product_has_new_N_subst:
                                    has_nucleophilic_coupling = True
                                    coupling_depth = depth
                                    if "N-functionalization" not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append("N-functionalization")
                                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                                    if checker.check_fg("Secondary amine", product) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                                    if checker.check_fg("Tertiary amine", product) and "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                except Exception:
                    pass

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If we found pre-activated heterocycles but no explicit activation reaction,
    # the sequence is correct if there's a coupling reaction.
    if pre_activated_heterocycles and activation_depth == -1 and coupling_depth != -1:
        reaction_sequence_correct = True
    # Otherwise, check if activation occurs before coupling in the synthesis
    # (higher depth means earlier in forward synthesis)
    elif activation_depth != -1 and coupling_depth != -1:
        if activation_depth > coupling_depth:
            reaction_sequence_correct = True

    strategy_present = (
        has_heterocycle
        and has_electrophilic_activation
        and has_nucleophilic_coupling
        and reaction_sequence_correct
    )

    # Record structural constraints if met
    if has_heterocycle and has_electrophilic_activation and has_nucleophilic_coupling:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain a heterocycle of interest, an electrophilic activation event, and a nucleophilic coupling event.",
                "targets": [
                    "HETEROCYCLE_OF_INTEREST",
                    "ELECTROPHILIC_ACTIVATION",
                    "NUCLEOPHILIC_COUPLING"
                ]
            }
        })
    
    if reaction_sequence_correct:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "description": "The electrophilic activation event (either a reaction or the use of a pre-activated starting material) must occur before the nucleophilic coupling reaction.",
                "first_event": "ELECTROPHILIC_ACTIVATION",
                "second_event": "NUCLEOPHILIC_COUPLING"
            }
        })

    return strategy_present, findings_json
