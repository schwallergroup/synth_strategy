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
    Detects if the synthetic route involves the formation of a stereocenter
    through reduction of a carbonyl group to an alcohol.
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

    # Track if we found the pattern
    found_reduction = False

    def check_stereocenter_formation(reactants_smiles, product_smiles, rsmi):
        nonlocal found_reduction, findings_json

        has_carbonyl = False
        for r in reactants_smiles:
            if checker.check_fg("Ketone", r):
                has_carbonyl = True
                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
            if checker.check_fg("Aldehyde", r):
                has_carbonyl = True
                if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")

        has_alcohol = False
        if checker.check_fg("Primary alcohol", product_smiles):
            has_alcohol = True
            if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
        if checker.check_fg("Secondary alcohol", product_smiles):
            has_alcohol = True
            if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
        if checker.check_fg("Tertiary alcohol", product_smiles):
            has_alcohol = True
            if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")

        if has_carbonyl and has_alcohol:
            try:
                reactants_mol = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if not all(reactants_mol) or not product_mol:
                    return

                potential_stereo = Chem.FindPotentialStereo(product_mol)

                carbonyl_carbons_with_map = []
                for mol_idx, mol in enumerate(reactants_mol):
                    if mol:
                        for atom in mol.GetAtoms():
                            if atom.GetSymbol() == "C" and atom.GetAtomMapNum() > 0:
                                is_carbonyl = False
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.GetSymbol() == "O":
                                        bond = mol.GetBondBetweenAtoms(
                                            atom.GetIdx(), neighbor.GetIdx()
                                        )
                                        if (
                                            bond
                                            and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                                        ):
                                            is_carbonyl = True
                                            break

                                if is_carbonyl:
                                    carbonyl_carbons_with_map.append(
                                        (mol_idx, atom.GetIdx(), atom.GetAtomMapNum())
                                    )

                alcohol_carbons_with_map = []
                if product_mol:
                    for atom in product_mol.GetAtoms():
                        if atom.GetSymbol() == "C" and atom.GetAtomMapNum() > 0:
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetSymbol() == "O":
                                    if neighbor.GetTotalNumHs() > 0:
                                        alcohol_carbons_with_map.append(
                                            (atom.GetIdx(), atom.GetAtomMapNum())
                                        )
                                        break

                for _, _, carbonyl_map in carbonyl_carbons_with_map:
                    for alcohol_idx, alcohol_map in alcohol_carbons_with_map:
                        if carbonyl_map == alcohol_map:
                            carbon = product_mol.GetAtomWithIdx(alcohol_idx)

                            has_chirality = (
                                carbon.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED
                            )

                            is_potential_stereo = False
                            for stereo in potential_stereo:
                                if stereo.centeredOn == alcohol_idx:
                                    is_potential_stereo = True
                                    break

                            was_stereocenter = False
                            for mol_idx, c_idx, c_map in carbonyl_carbons_with_map:
                                if c_map == alcohol_map:
                                    reactant_carbon = reactants_mol[mol_idx].GetAtomWithIdx(
                                        c_idx
                                    )
                                    if (
                                        reactant_carbon.GetChiralTag()
                                        != Chem.rdchem.ChiralType.CHI_UNSPECIFIED
                                    ):
                                        was_stereocenter = True
                                        break

                            if not was_stereocenter and (
                                has_chirality or is_potential_stereo
                            ):
                                found_reduction = True
                                # Record the structural constraint
                                if {"type": "co-occurrence", "details": {"description": "A single reaction step must be a recognized carbonyl reduction, have a carbonyl group in the reactants, and an alcohol group in the product, where the transformation creates a new stereocenter."}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "A single reaction step must be a recognized carbonyl reduction, have a carbonyl group in the reactants, and an alcohol group in the product, where the transformation creates a new stereocenter."}})
                                return

            except Exception as e:
                pass

    def dfs_traverse(node, depth=0):
        nonlocal found_reduction, findings_json

        if found_reduction:
            return  # Early exit if we already found what we're looking for

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            is_reduction = False
            reduction_types = [
                "Reduction of aldehydes and ketones to alcohols",
                "Reduction of ketone to secondary alcohol",
                "Grignard from ketone to alcohol",
                "Grignard from aldehyde to alcohol"
            ]

            for r_name in reduction_types:
                if checker.check_reaction(r_name, rsmi):
                    is_reduction = True
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break

            if is_reduction:
                check_stereocenter_formation(reactants_smiles, product_smiles, rsmi)

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return found_reduction, findings_json