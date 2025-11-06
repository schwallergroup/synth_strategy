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


# Refactoring for Enumeration: Isolate the list of heterocycles.
HETEROCYCLES_OF_INTEREST = [
    "pyridine",
    "pyrrole",
    "furan",
    "thiophene",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrazole",
    "isoxazole",
    "isothiazole",
    "triazole",
    "tetrazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "indole",
    "benzofuran",
    "benzothiophene",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "quinoline",
    "isoquinoline",
    "purine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy for assembling molecules containing specific heterocyclic scaffolds. A route is identified as following this strategy if it meets three criteria: 1) It uses at least one key C-C or C-N bond-forming reaction (specifically, Suzuki or amide coupling). 2) It involves heterocycles from the HETEROCYCLES_OF_INTEREST list, either by forming one or using one as a building block, or includes a nitro group reduction. 3) It is convergent, defined as combining at least three unique fragments throughout the synthesis.
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
    has_suzuki = False
    has_amide_formation = False
    has_nitro_reduction = False
    heterocycle_types = set()
    unique_fragments = set()
    heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki, has_amide_formation, has_nitro_reduction, heterocycle_types, unique_fragments, heterocycle_formation, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for heterocycles in molecules
            for heterocycle in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(heterocycle, mol_smiles):
                    heterocycle_types.add(heterocycle)
                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                    # print(f"Found heterocycle: {heterocycle} in {mol_smiles}")

        elif node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Add reactants to unique fragments set (using canonical SMILES)
            for reactant in reactants_smiles:
                try:
                    canonical_smiles = Chem.CanonSmiles(reactant)
                    unique_fragments.add(canonical_smiles)
                except:
                    unique_fragments.add(reactant)  # Fallback if canonicalization fails

            # Check for Suzuki coupling (more comprehensive)
            suzuki_reactions = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki",  # Generic Suzuki check
            ]

            for rxn_type in suzuki_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    has_suzuki = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    # print(f"Found Suzuki coupling at depth {depth}: {rxn_type}")
                    break

            # Check for amide formation (expanded)
            amide_reactions = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Schotten-Baumann_amide",
            ]

            for rxn_type in amide_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    has_amide_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    # print(f"Found amide formation at depth {depth}: {rxn_type}")
                    break

            # Check for nitro reduction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                has_nitro_reduction = True
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                # print(f"Found nitro reduction at depth {depth}")

            # Check for heterocycle formation
            reactant_heterocycles = set()
            product_heterocycles = set()

            for heterocycle in HETEROCYCLES_OF_INTEREST:
                for reactant in reactants_smiles:
                    if checker.check_ring(heterocycle, reactant):
                        reactant_heterocycles.add(heterocycle)

                if checker.check_ring(heterocycle, product_smiles):
                    product_heterocycles.add(heterocycle)
                    heterocycle_types.add(heterocycle)
                    # print(f"Found heterocycle: {heterocycle} in product {product_smiles}")

            # Check if a new heterocycle was formed
            new_heterocycles = product_heterocycles - reactant_heterocycles
            if new_heterocycles:
                heterocycle_formation = True
                # Add 'ring_formation' to named_reactions if a new heterocycle was formed
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                for hc in new_heterocycles:
                    if hc not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(hc)
                # print(f"Heterocycle formation detected at depth {depth}: {new_heterocycles}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if the route matches the convergent heterocycle assembly strategy
    is_convergent = (
        (has_suzuki or has_amide_formation)  # At least one coupling reaction type
        and (
            has_nitro_reduction or heterocycle_formation or len(heterocycle_types) >= 1
        )  # Either nitro reduction, heterocycle formation, or presence of heterocycles
        and len(unique_fragments) >= 3  # At least 3 unique fragments (convergent)
    )

    # Populate structural constraints based on the final flags
    if has_suzuki or has_amide_formation:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "suzuki_coupling_group",
                    "amide_formation_group"
                ],
                "min_required": 1,
                "description": "The route must contain at least one reaction from the Suzuki coupling group or the amide formation group."
            }
        })

    if has_nitro_reduction or heterocycle_formation or len(heterocycle_types) >= 1:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Reduction of nitro groups to amines",
                    "ring_formation",
                    "heterocycle_presence"
                ],
                "min_required": 1,
                "description": "The route must contain either a nitro reduction, a de novo formation of a specified heterocycle, or the presence of at least one specified heterocycle type as a building block."
            }
        })

    if len(unique_fragments) >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_fragments",
                "operator": ">=",
                "value": 3,
                "description": "The route must combine at least 3 unique fragments, indicating a convergent strategy."
            }
        })

    # print(f"Convergent heterocycle assembly analysis:")
    # print(f"- Suzuki coupling: {has_suzuki}")
    # print(f"- Amide formation: {has_amide_formation}")
    # print(f"- Nitro reduction: {has_nitro_reduction}")
    # print(f"- Heterocycle formation: {heterocycle_formation}")
    # print(f"- Heterocycle types: {len(heterocycle_types)} {heterocycle_types}")
    # print(f"- Unique fragments: {len(unique_fragments)}")
    # print(f"- Is convergent heterocycle assembly: {is_convergent}")

    return is_convergent, findings_json
