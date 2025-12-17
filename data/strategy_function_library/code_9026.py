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


# Refactoring for Enumeration: Isolate lists of interest
HETEROCYCLIC_RINGS_OF_INTEREST = [
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole", "indole",
    "quinoline", "isoquinoline", "purine", "benzimidazole", "piperidine",
    "piperazine", "morpholine", "thiomorpholine", "pyrrolidine",
]

N_FUNCTIONALIZATION_REACTIONS = [
    "N-methylation", "Methylation with MeI_primary", "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary", "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde", "DMS Amine methylation",
    "N-alkylation", "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage N-functionalization (e.g., alkylation, arylation) of specific heterocyclic scaffolds. The strategy is confirmed by identifying a relevant reaction from `N_FUNCTIONALIZATION_REACTIONS` and verifying that a reactant containing a heterocyclic N-H group (from `HETEROCYCLIC_RINGS_OF_INTEREST`) is consumed.
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

    n_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected, findings_json

        if n_alkylation_detected:
            return

        if node["type"] == "reaction" and depth <= 3:
            # Record positional constraint if met
            if depth <= 3:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "N-functionalization of heterocycle",
                        "position": "late_stage",
                        "max_depth": 3
                    }
                })

            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a relevant N-functionalization reaction
                is_n_functionalization = False
                for rxn_type in N_FUNCTIONALIZATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_n_functionalization = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_n_functionalization:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if the product contains a heterocycle of interest
                    product_has_heterocycle = False
                    for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                        if checker.check_ring(ring, product):
                            product_has_heterocycle = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            break
                    
                    if product_has_heterocycle:
                        # Check if any reactant has a heterocycle with N-H
                        reactant_has_target_heterocycle_nh = False
                        for reactant in reactants:
                            reactant_has_target_ring = False
                            for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                                if checker.check_ring(ring, reactant):
                                    reactant_has_target_ring = True
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                                    break
                            
                            if reactant_has_target_ring:
                                try:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol:
                                        # Check specifically for a heterocyclic N-H
                                        has_heterocyclic_nh = False
                                        for atom in reactant_mol.GetAtoms():
                                            if (
                                                atom.GetSymbol() == "N"
                                                and atom.IsInRing()
                                                and atom.GetTotalNumHs() > 0
                                            ):
                                                has_heterocyclic_nh = True
                                                break
                                        
                                        if has_heterocyclic_nh:
                                            reactant_has_target_heterocycle_nh = True
                                            if "heterocyclic N-H" not in findings_json["atomic_checks"]["functional_groups"]:
                                                findings_json["atomic_checks"]["functional_groups"].append("heterocyclic N-H")
                                            break # Found a reactant with N-H heterocycle
                                except Exception:
                                    # Error processing reactant SMILES, continue
                                    pass
                        
                        if is_n_functionalization and product_has_heterocycle and reactant_has_target_heterocycle_nh:
                            n_alkylation_detected = True
                            # Record co-occurrence constraint if all conditions met
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "N-functionalization reaction",
                                        "reactant with target heterocycle and N-H",
                                        "product with target heterocycle"
                                    ],
                                    "scope": "reaction_step"
                                }
                            })
                            return

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return n_alkylation_detected, findings_json
