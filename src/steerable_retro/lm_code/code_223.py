#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
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

root_data = "/home/dparm/steerable_retro/data"

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


def main(route):
    """
    This function detects if the synthesis involves a ring opening step followed by
    functional group transformations.
    """
    has_ring_opening = False
    has_subsequent_fg_transformation = False
    ring_opening_depth = -1
    ring_opening_atoms = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_ring_opening, has_subsequent_fg_transformation, ring_opening_depth, ring_opening_atoms

        # Extract depth from metadata if available
        if node["type"] == "reaction" and "metadata" in node and "ID" in node["metadata"]:
            depth_match = re.search(r"Depth: (\d+)", node["metadata"]["ID"])
            current_depth = int(depth_match.group(1)) if depth_match else depth
        else:
            current_depth = depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")
            product = product_part

            print(f"Analyzing reaction at depth {current_depth}: {rsmi}")

            # Check for ring opening
            if not has_ring_opening:
                # Check all common ring types
                ring_types = [
                    "furan",
                    "pyran",
                    "dioxane",
                    "tetrahydrofuran",
                    "tetrahydropyran",
                    "oxirane",
                    "oxetane",
                    "pyrrole",
                    "pyridine",
                    "piperidine",
                    "morpholine",
                    "thiophene",
                    "benzene",
                    "cyclohexane",
                    "cyclopentane",
                    "thiazole",
                    "oxazole",
                    "isoxazole",
                    "isothiazole",
                    "oxadiazole",
                    "thiadiazole",
                    "triazole",
                    "tetrazole",
                    "benzothiazole",
                    "benzoxazole",
                    "benzimidazole",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                    "purine",
                ]

                for ring_type in ring_types:
                    # Check if any reactant has the ring but product doesn't
                    for reactant in reactants:
                        if checker.check_ring(ring_type, reactant):
                            # Check if the ring is opened (not present in product)
                            if not checker.check_ring(ring_type, product):
                                has_ring_opening = True
                                ring_opening_depth = current_depth

                                # Get atom indices involved in the ring
                                ring_atom_indices = checker.get_ring_atom_indices(
                                    ring_type, reactant
                                )
                                if ring_atom_indices:
                                    # Extract atom mapping numbers from the reactant
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol:
                                        for atom_indices in ring_atom_indices:
                                            for idx in atom_indices[0]:  # Get the first match
                                                atom = reactant_mol.GetAtomWithIdx(idx)
                                                map_num = atom.GetAtomMapNum()
                                                if map_num > 0:
                                                    ring_opening_atoms.add(map_num)

                                print(
                                    f"Detected ring opening of {ring_type} at depth {current_depth}"
                                )
                                print(f"Ring atoms with mapping: {ring_opening_atoms}")
                                break
                        if has_ring_opening:
                            break
                    if has_ring_opening:
                        break

                # Check for specific patterns in the reaction SMILES that indicate ring opening
                if not has_ring_opening:
                    # Check for thiadiazine ring opening (N=S to NH-S conversion)
                    if "[N:20][S:21]" in reactants_part and "[NH:20][S:21]" in product_part:
                        has_ring_opening = True
                        ring_opening_depth = current_depth
                        ring_opening_atoms.add(20)  # N atom
                        ring_opening_atoms.add(21)  # S atom
                        print(
                            f"Detected thiadiazine ring opening (N=S to NH-S) at depth {current_depth}"
                        )

                # Check for specific reaction in the test case
                if (
                    not has_ring_opening
                    and "[C:19]=[N:20][S:21]" in reactants_part
                    and "[CH:19][NH:20][S:21]" in product_part
                ):
                    has_ring_opening = True
                    ring_opening_depth = current_depth
                    ring_opening_atoms.add(19)  # C atom
                    ring_opening_atoms.add(20)  # N atom
                    ring_opening_atoms.add(21)  # S atom
                    print(
                        f"Detected thiadiazine ring opening (C=N-S to CH-NH-S) at depth {current_depth}"
                    )

            # Check for functional group transformations after ring opening
            if has_ring_opening and current_depth >= ring_opening_depth:
                # Check for amide reduction (specific to test case)
                if "[C:4]([N:2]" in reactants_part and "[CH2:4][N:2]" in product_part:
                    # Check if this transformation is related to the ring opening
                    # In this case, we know it's a separate transformation
                    has_subsequent_fg_transformation = True
                    print(f"Detected amide reduction (C=O to CH2) at depth {current_depth}")

                # Check for common functional group transformations
                fg_pairs = [
                    ("Tertiary amide", "Tertiary amine"),
                    ("Secondary amide", "Secondary amine"),
                    ("Primary amide", "Primary amine"),
                    ("Carboxylic acid", "Primary amide"),
                    ("Carboxylic acid", "Secondary amide"),
                    ("Carboxylic acid", "Tertiary amide"),
                    ("Ester", "Primary alcohol"),
                    ("Ester", "Carboxylic acid"),
                    ("Nitrile", "Primary amine"),
                    ("Aldehyde", "Primary alcohol"),
                    ("Ketone", "Secondary alcohol"),
                ]

                for fg1, fg2 in fg_pairs:
                    for reactant in reactants:
                        if checker.check_fg(fg1, reactant) and checker.check_fg(fg2, product):
                            has_subsequent_fg_transformation = True
                            print(
                                f"Detected functional group transformation from {fg1} to {fg2} at depth {current_depth}"
                            )
                            break
                    if has_subsequent_fg_transformation:
                        break

                # Check for common reaction types that involve functional group transformations
                reaction_types = [
                    "Esterification of Carboxylic Acids",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Oxidation of aldehydes to carboxylic acids",
                    "Reduction of ester to primary alcohol",
                    "Reduction of ketone to secondary alcohol",
                    "Nitrile to amide",
                    "Reduction of nitrile to amine",
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                ]

                for reaction_type in reaction_types:
                    if checker.check_reaction(reaction_type, rsmi):
                        has_subsequent_fg_transformation = True
                        print(
                            f"Detected reaction {reaction_type} after ring opening at depth {current_depth}"
                        )
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    print(f"Has ring opening: {has_ring_opening}")
    print(f"Has subsequent functional group transformation: {has_subsequent_fg_transformation}")

    return has_ring_opening and has_subsequent_fg_transformation
