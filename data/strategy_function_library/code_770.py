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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects the use of N-bromosuccinimide (NBS) as a brominating agent. The strategy is confirmed if NBS is a reactant and: 1) the reaction is a known bromination type, 2) succinimide is a product, or 3) the number of bromine atoms increases in the product.
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

    nbs_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nbs_detected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Multiple representations of NBS
            nbs_smiles_variants = ["O=C1CCC(=O)N1Br", "BrN1C(=O)CCC1=O", "O=C1N(Br)C(=O)CCC1"]
            succinimide_smiles_variants = ["O=C1CCC(=O)N1", "N1C(=O)CCC1=O"]

            # Check if any reactant is NBS using substructure matching
            nbs_present = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant.strip())
                if reactant_mol:
                    for nbs_variant in nbs_smiles_variants:
                        nbs_mol = Chem.MolFromSmiles(nbs_variant)
                        if nbs_mol and (
                            Chem.MolToSmiles(reactant_mol) == Chem.MolToSmiles(nbs_mol)
                            or reactant_mol.HasSubstructMatch(nbs_mol)
                        ):
                            nbs_present = True
                            if "N-bromosuccinimide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("N-bromosuccinimide")
                            print(f"Found NBS in reactants: {reactant}")
                            break
                if nbs_present:
                    break

            # Check for Wohl-Ziegler and other bromination reactions
            bromination_reaction_types = [
                "Wohl-Ziegler bromination benzyl primary",
                "Wohl-Ziegler bromination benzyl secondary",
                "Wohl-Ziegler bromination benzyl tertiary",
                "Wohl-Ziegler bromination allyl primary",
                "Wohl-Ziegler bromination allyl secondary",
                "Wohl-Ziegler bromination allyl tertiary",
                "Wohl-Ziegler bromination carbonyl primary",
                "Wohl-Ziegler bromination carbonyl secondary",
                "Wohl-Ziegler bromination carbonyl tertiary",
                "Aromatic bromination",
                "Bromination",
            ]

            is_bromination_reaction = False
            detected_bromination_type = None
            for rxn_type in bromination_reaction_types:
                if checker.check_reaction(rxn_type, rsmi):
                    is_bromination_reaction = True
                    detected_bromination_type = rxn_type
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            # Check for succinimide as a byproduct
            succinimide_present = False
            for product_part in product.split("."):
                product_part_mol = Chem.MolFromSmiles(product_part.strip())
                if product_part_mol:
                    for succinimide_variant in succinimide_smiles_variants:
                        succinimide_mol = Chem.MolFromSmiles(succinimide_variant)
                        if succinimide_mol and (
                            Chem.MolToSmiles(product_part_mol) == Chem.MolToSmiles(succinimide_mol)
                            or product_part_mol.HasSubstructMatch(succinimide_mol)
                        ):
                            succinimide_present = True
                            if "succinimide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("succinimide")
                            print(f"Found succinimide in products: {product_part}")
                            break
                if succinimide_present:
                    break

            if nbs_present and is_bromination_reaction:
                print(f"Found NBS bromination reaction: {rsmi}")
                nbs_detected = True
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "N-bromosuccinimide",
                            "named_bromination_reaction"
                        ]
                    }
                })
            elif nbs_present and succinimide_present:
                print(f"Found NBS bromination with succinimide byproduct: {rsmi}")
                nbs_detected = True
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "N-bromosuccinimide",
                            "succinimide"
                        ]
                    }
                })
            # If we didn't detect a specific bromination reaction but NBS is present,
            # check for bromination by comparing functional groups
            elif nbs_present:
                print(f"Found NBS reagent, checking for bromination: {rsmi}")

                # Get non-NBS reactants
                non_nbs_reactants = []
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant.strip())
                    if reactant_mol:
                        is_nbs = False
                        for nbs_variant in nbs_smiles_variants:
                            nbs_mol = Chem.MolFromSmiles(nbs_variant)
                            if nbs_mol and (
                                Chem.MolToSmiles(reactant_mol) == Chem.MolToSmiles(nbs_mol)
                                or reactant_mol.HasSubstructMatch(nbs_mol)
                            ):
                                is_nbs = True
                                break
                        if not is_nbs:
                            non_nbs_reactants.append(reactant)

                if non_nbs_reactants:
                    # Count bromine atoms properly using RDKit
                    def count_bromine_atoms(smiles):
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 35)
                        return 0

                    reactant_br_count = sum(count_bromine_atoms(r) for r in non_nbs_reactants)
                    product_br_count = count_bromine_atoms(product)

                    if product_br_count > reactant_br_count:
                        print(f"Detected bromination with NBS: {rsmi}")
                        nbs_detected = True
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "bromine_atoms",
                                "operator": ">",
                                "value": "reactant_non_nbs_count"
                            }
                        })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return nbs_detected, findings_json