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


CYCLIZATION_REACTIONS = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "Benzothiazole formation from aldehyde", "Benzothiazole formation from acyl halide",
    "Benzothiazole formation from ester/carboxylic acid", "Benzoxazole formation from aldehyde",
    "Benzoxazole formation from acyl halide", "Benzoxazole formation from ester/carboxylic acid",
    "Benzoxazole formation (intramolecular)", "Benzimidazole formation from aldehyde",
    "Benzimidazole formation from acyl halide", "Benzimidazole formation from ester/carboxylic acid",
    "Intramolecular amination (heterocycle formation)",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)", "Diels-Alder",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition", "Pyrazole formation",
    "Michael-induced ring closure from hydrazone", "Michael-induced ring closure from diazoalkane",
    "[3+2]-cycloaddition of hydrazone and alkyne", "[3+2]-cycloaddition of hydrazone and alkene",
    "[3+2]-cycloaddition of diazoalkane and alkyne", "[3+2]-cycloaddition of diazoalkane and alkene",
    "Pictet-Spengler",
]

HETEROCYCLIC_RINGS_OF_INTEREST = [
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "benzimidazole", "benzoxazole", "benzothiazole", "indole", "quinoline",
    "isoquinoline", "furan", "thiophene", "triazole", "tetrazole", "isoxazole",
    "isothiazole", "oxadiazole", "thiadiazole",
]

AROMATIC_RINGS_OF_INTEREST = [
    "benzene", "pyridine", "furan", "thiophene", "pyrrole", "imidazole",
    "oxazole", "thiazole", "indole", "quinoline", "isoquinoline", "naphthalene",
]

FUNCTIONALIZATION_REACTIONS = [
    "Aromatic bromination", "Aromatic chlorination", "Aromatic fluorination",
    "Aromatic iodination", "Aromatic nitration with HNO3", "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt", "Aromatic nitration with alkyl NO2",
    "Friedel-Crafts acylation", "Friedel-Crafts alkylation",
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Heck terminal vinyl", "Heck reaction with vinyl ester and amine",
    "Negishi coupling", "Stille reaction_aryl", "Hiyama-Denmark Coupling",
    "Kumada cross-coupling", "Aryllithium cross-coupling",
    "Directed ortho metalation of arenes", "Minisci (para)", "Minisci (ortho)",
    "Catellani reaction ortho", "Catellani reaction para",
]

AROMATIC_FGS_OF_INTEREST = [
    "Aromatic halide", "Nitro group", "Phenol", "Aniline", "Triflate",
    "Mesylate", "Tosylate", "Boronic acid", "Boronic ester", "Nitrile",
    "Carboxylic acid", "Ester", "Aldehyde", "Ketone", "Primary amide",
    "Secondary amide", "Tertiary amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a 'build-and-cyclize' strategy, defined by early-stage aromatic functionalization (depth > 2) 
    followed by late-stage cyclization (depth <= 2). Aromatic functionalization is identified by specific 
    named reactions (see FUNCTIONALIZATION_REACTIONS) or the addition of key functional groups (see 
    AROMATIC_FGS_OF_INTEREST) to an aromatic ring (see AROMATIC_RINGS_OF_INTEREST). Cyclization is 
    identified by an increase in ring count, specific named reactions (see CYCLIZATION_REACTIONS), or the 
    formation of a new heterocycle (see HETEROCYCLIC_RINGS_OF_INTEREST).
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

    # Track key transformations and their depths
    aromatic_functionalizations = []
    cyclizations = []

    def dfs_traverse(node, current_depth=0):
        nonlocal aromatic_functionalizations, cyclizations, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            depth = node.get("metadata", {}).get("depth", current_depth)

            if depth == current_depth and "ID" in node.get("metadata", {}):
                depth_match = None
                id_value = node["metadata"]["ID"]
                if isinstance(id_value, str):
                    depth_match = id_value.split("Depth: ")
                    if len(depth_match) > 1:
                        try:
                            depth = int(depth_match[1].split()[0])
                        except (ValueError, IndexError):
                            pass

            try:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                ]

                if product_mol and reactant_mols:
                    product_ring_count = len(Chem.GetSSSR(product_mol))
                    reactant_ring_counts = [len(Chem.GetSSSR(r)) for r in reactant_mols]
                    max_reactant_ring_count = (
                        max(reactant_ring_counts) if reactant_ring_counts else 0
                    )

                    is_cyclization_reaction = False
                    for rxn in CYCLIZATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_cyclization_reaction = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)

                    ring_formed = product_ring_count > max_reactant_ring_count
                    if ring_formed:
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                    new_ring_formed = False
                    for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                        if checker.check_ring(ring, product_smiles):
                            if not any(checker.check_ring(ring, r) for r in reactants_smiles):
                                new_ring_formed = True
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)

                    if ring_formed or is_cyclization_reaction or new_ring_formed:
                        cyclizations.append(depth)

                has_aromatic_reactant = False
                for r in reactants_smiles:
                    for ring in AROMATIC_RINGS_OF_INTEREST:
                        if checker.check_ring(ring, r):
                            has_aromatic_reactant = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            break
                    if has_aromatic_reactant: break

                is_functionalization = False
                for rxn in FUNCTIONALIZATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_functionalization = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                has_aromatic_fg_in_product = False
                for fg in AROMATIC_FGS_OF_INTEREST:
                    if checker.check_fg(fg, product_smiles):
                        has_aromatic_fg_in_product = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

                new_fg_added = False
                if has_aromatic_fg_in_product:
                    for fg in AROMATIC_FGS_OF_INTEREST:
                        if checker.check_fg(fg, product_smiles):
                            if not any(checker.check_fg(fg, r) for r in reactants_smiles):
                                new_fg_added = True

                if has_aromatic_reactant and (is_functionalization or new_fg_added):
                    aromatic_functionalizations.append(depth)

            except Exception:
                pass

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (chemicals)
                next_depth = current_depth
            else:
                # If current node is a chemical, depth increases for its children (reactions)
                next_depth = current_depth + 1
            dfs_traverse(child, next_depth)

    dfs_traverse(route)

    has_early_functionalization = any(depth > 2 for depth in aromatic_functionalizations)
    has_late_cyclization = any(depth <= 2 for depth in cyclizations)

    result = has_early_functionalization and has_late_cyclization

    if result:
        # Add the structural constraint if the overall condition is met
        structural_constraint_obj = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "early_stage_aromatic_functionalization",
                    "late_stage_cyclization"
                ]
            }
        }
        if structural_constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(structural_constraint_obj)

    return result, findings_json
