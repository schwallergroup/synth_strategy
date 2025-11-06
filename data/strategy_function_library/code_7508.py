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


# Refactoring for Enumeration: Isolate lists of chemical entities
DIAZINES_OF_INTEREST = ["pyrimidine", "pyridazine", "pyrazine"]

TRIAZOLE_FORMATION_REACTIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "{Huisgen_Cu-catalyzed_1,4-subst}",
    "{Huisgen_Ru-catalyzed_1,5_subst}",
    "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}",
    "Azide-nitrile click cycloaddition to triazole",
    "{1,2,4-triazole_acetohydrazide}",
    "{1,2,4-triazole_carboxylic-acid/ester}",
]

N_ALKYLATION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Alkylation of amines",
    "N-methylation",
    "{Mitsunobu_tetrazole_1}",
    "{Mitsunobu_tetrazole_2}",
    "{Mitsunobu_tetrazole_3}",
    "{Mitsunobu_tetrazole_4}",
    "Mitsunobu_imide",
]

FRAGMENT_COUPLING_REACTIONS = [
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "{Buchwald-Hartwig}",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "{heteroaromatic_nuc_sub}",
    "Ullmann-Goldberg Substitution amine",
    "Goldberg coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a linear build-up strategy for constructing
    a complex heterocyclic scaffold (triazolopyrimidine).
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

    # Track the build-up steps
    heterocycle_formation = False
    n_alkylation = False
    fragment_coupling = False
    final_scaffold_present = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation, n_alkylation, fragment_coupling, final_scaffold_present, findings_json

        # Check if the final product contains the target scaffolds
        if node["type"] == "mol" and depth == 0:
            mol_smiles = node["smiles"]
            has_triazole = checker.check_ring("triazole", mol_smiles)
            if has_triazole:
                findings_json["atomic_checks"]["ring_systems"].append("triazole")

            has_diazine = False
            for d in DIAZINES_OF_INTEREST:
                if checker.check_ring(d, mol_smiles):
                    has_diazine = True
                    findings_json["atomic_checks"]["ring_systems"].append(d)

            if has_triazole and has_diazine:
                final_scaffold_present = True
                # This corresponds to the 'positional' structural constraint
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "co-occurrence of triazole and a diazine ring (pyrimidine, pyridazine, or pyrazine)",
                        "position": "last_stage"
                    }
                })

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                try:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check for triazole formation
                    product_has_triazole = checker.check_ring("triazole", product_smiles)
                    reactants_have_triazole = any(
                        checker.check_ring("triazole", r) for r in reactants_smiles
                    )
                    if product_has_triazole and not reactants_have_triazole:
                        heterocycle_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                    # Check for N-alkylation
                    for rxn in N_ALKYLATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            n_alkylation = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                    # Check for fragment coupling
                    for rxn in FRAGMENT_COUPLING_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            fragment_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                except Exception:
                    # Silently ignore errors in reaction processing
                    pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or other types that should increase depth
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = heterocycle_formation and n_alkylation and fragment_coupling and final_scaffold_present

    if result:
        # This corresponds to the 'co-occurrence' structural constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "triazole_ring_formation",
                    "N-alkylation_reaction",
                    "fragment_coupling_reaction",
                    "final_product_scaffold"
                ]
            }
        })

    return result, findings_json
