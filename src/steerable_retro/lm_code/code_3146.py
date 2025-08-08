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
    This function detects if the synthetic route involves formation of aromatic C-N bonds.
    """
    has_aromatic_cn_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aromatic_cn_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction: {rsmi}")

            # Check for specific N-arylation reactions
            n_arylation_reactions = [
                "Buchwald-Hartwig",
                "N-arylation",
                "Ullmann-Goldberg",
                "Goldberg coupling",
                "N-arylation_heterocycles",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "Goldberg coupling aryl amine-aryl chloride",
                "Goldberg coupling aryl amide-aryl chloride",
                "Ullmann-Goldberg Substitution amine",
                "Catellani reaction ortho",
                "Catellani reaction para",
                "Minisci (para)",
                "Minisci (ortho)",
                "Minisci (para-cyanide)",
                "Minisci (ortho-cyanide)",
                "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                "Intramolecular amination (heterocycle formation)",
            ]

            for reaction_type in n_arylation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected {reaction_type} reaction: {rsmi}")
                    has_aromatic_cn_formation = True
                    return

            # Check for heterocycle formation reactions that create aromatic C-N bonds
            heterocycle_formation_reactions = [
                "Formation of NOS Heterocycles",
                "benzimidazole_derivatives_carboxylic-acid/ester",
                "benzimidazole_derivatives_aldehyde",
                "benzothiazole",
                "benzoxazole_arom-aldehyde",
                "benzoxazole_carboxylic-acid",
                "thiazole",
                "Niementowski_quinazoline",
                "tetrazole_terminal",
                "tetrazole_connect_regioisomere_1",
                "tetrazole_connect_regioisomere_2",
                "1,2,4-triazole_acetohydrazide",
                "1,2,4-triazole_carboxylic-acid/ester",
                "3-nitrile-pyridine",
                "pyrazole",
                "Paal-Knorr pyrrole",
                "triaryl-imidazole",
                "Fischer indole",
                "Friedlaender chinoline",
                "indole",
                "oxadiazole",
                "imidazole",
            ]

            for reaction_type in heterocycle_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected heterocycle formation reaction: {reaction_type}")
                    has_aromatic_cn_formation = True
                    return

            # If no specific reaction type matched, check for aromatic C-N bond formation
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aromatic C-N bond formation using atom mapping
            product_mol = Chem.MolFromSmiles(product_smiles)
            if not product_mol:
                for child in node.get("children", []):
                    dfs_traverse(child, depth + 1)
                return

            # Check for nitrogen-containing heterocycles in product but not in reactants
            nitrogen_heterocycles = [
                "pyridine",
                "pyrazole",
                "imidazole",
                "oxazole",
                "thiazole",
                "pyrimidine",
                "pyrazine",
                "pyridazine",
                "triazole",
                "tetrazole",
                "indole",
                "quinoline",
                "isoquinoline",
                "purine",
                "carbazole",
                "acridine",
                "benzimidazole",
                "benzoxazole",
                "benzothiazole",
            ]

            for ring in nitrogen_heterocycles:
                if checker.check_ring(ring, product_smiles):
                    ring_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_ring(ring, reactant):
                            ring_in_reactants = True
                            break

                    if not ring_in_reactants:
                        print(f"Detected formation of nitrogen heterocycle: {ring}")
                        has_aromatic_cn_formation = True
                        return

            # Look for aromatic carbon - nitrogen bonds using multiple patterns
            aromatic_cn_patterns = [
                Chem.MolFromSmarts("c-[#7]"),  # Single bond
                Chem.MolFromSmarts("c:[#7]"),  # Aromatic bond
                Chem.MolFromSmarts("c~[#7]"),  # Any bond
            ]

            for pattern in aromatic_cn_patterns:
                if not pattern:
                    continue

                product_matches = product_mol.GetSubstructMatches(pattern)

                # Check each bond in product to see if it's new
                for match in product_matches:
                    c_atom_idx, n_atom_idx = match

                    # Get atom maps if available
                    c_atom_map = (
                        product_mol.GetAtomWithIdx(c_atom_idx).GetProp("molAtomMapNumber")
                        if product_mol.GetAtomWithIdx(c_atom_idx).HasProp("molAtomMapNumber")
                        else None
                    )
                    n_atom_map = (
                        product_mol.GetAtomWithIdx(n_atom_idx).GetProp("molAtomMapNumber")
                        if product_mol.GetAtomWithIdx(n_atom_idx).HasProp("molAtomMapNumber")
                        else None
                    )

                    if c_atom_map and n_atom_map:
                        # Check if these mapped atoms are bonded in any reactant
                        bond_exists_in_reactants = False
                        for reactant in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if not reactant_mol:
                                continue

                            # Find atoms with these mapping numbers
                            c_reactant_idx = None
                            n_reactant_idx = None

                            for atom in reactant_mol.GetAtoms():
                                if atom.HasProp("molAtomMapNumber"):
                                    atom_map = atom.GetProp("molAtomMapNumber")
                                    if atom_map == c_atom_map:
                                        c_reactant_idx = atom.GetIdx()
                                    elif atom_map == n_atom_map:
                                        n_reactant_idx = atom.GetIdx()

                            # Check if these atoms exist and are bonded in this reactant
                            if c_reactant_idx is not None and n_reactant_idx is not None:
                                bond = reactant_mol.GetBondBetweenAtoms(
                                    c_reactant_idx, n_reactant_idx
                                )
                                if bond:
                                    bond_exists_in_reactants = True
                                    break

                        if not bond_exists_in_reactants:
                            print(
                                f"Detected new aromatic C-N bond formation between mapped atoms {c_atom_map} and {n_atom_map}"
                            )
                            has_aromatic_cn_formation = True
                            return

            # Check for aromatic amine formation
            aromatic_amine_fgs = ["Aniline", "Primary amine", "Secondary amine", "Tertiary amine"]
            for fg in aromatic_amine_fgs:
                if checker.check_fg(fg, product_smiles):
                    # Check if the aromatic amine is present in any reactant
                    fg_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_fg(fg, reactant):
                            fg_in_reactants = True
                            break

                    if not fg_in_reactants:
                        print(f"Detected formation of {fg}: {rsmi}")

                        # Verify it's an aromatic amine by checking if nitrogen is attached to aromatic carbon
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            aromatic_pattern = Chem.MolFromSmarts("c~[#7]")
                            if aromatic_pattern and product_mol.HasSubstructMatch(aromatic_pattern):
                                has_aromatic_cn_formation = True
                                return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {has_aromatic_cn_formation}")
    return has_aromatic_cn_formation
