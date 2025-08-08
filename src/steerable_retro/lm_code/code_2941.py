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
    This function detects if the synthesis follows a linear strategy
    where each step has only one non-reagent reactant.

    A linear synthesis strategy means that at each step, only one major building block
    is added to the growing molecule, while other reactants are considered reagents.
    """
    is_linear = True

    # Common reagent patterns
    common_reagents = [
        "H2O",
        "H2",
        "O2",
        "N2",
        "CO2",
        "CO",
        "HCl",
        "H2SO4",
        "HNO3",
        "NaOH",
        "KOH",
        "NH3",
        "CH3OH",
        "CH3CH2OH",
        "CH3COOH",
        "NaCl",
        "KCl",
        "NaBH4",
        "LiAlH4",
        "Pd/C",
        "SOCl2",
        "PCl3",
        "PCl5",
        "PBr3",
        "I2",
        "Br2",
        "Cl2",
        "NaH",
        "KH",
        "NaNH2",
        "LDA",
        "BuLi",
        "MeLi",
        "Et3N",
        "DIBAL",
        "DCC",
        "EDC",
        "HOBt",
        "DMAP",
        "TsCl",
        "MsCl",
        "Ac2O",
        "Boc2O",
        "TFA",
        "HBr",
        "HI",
        "H2O2",
        "KMnO4",
        "OsO4",
        "mCPBA",
        "NBS",
        "NCS",
        "NIS",
        "TBAF",
        "TMSCl",
        "TBSCl",
        "TBDPSCl",
        "MeI",
        "EtI",
        "BnBr",
        "AllylBr",
        "MgBr2",
        "ZnCl2",
        "CuI",
        "CuCl",
        "CuBr",
        "Pd(OAc)2",
        "Pd(PPh3)4",
        "Pd(dba)2",
        "Pt/C",
        "Raney Ni",
        "NaBH3CN",
        "NaCNBH3",
        "B2H6",
        "BH3",
        "BBr3",
        "BCl3",
        "AlCl3",
        "TiCl4",
        "SnCl2",
        "SnCl4",
        "FeCl3",
        "CrO3",
        "SeO2",
        "NaIO4",
        "LiOH",
        "K2CO3",
        "Na2CO3",
        "NaHCO3",
        "KHCO3",
        "MgSO4",
        "Na2SO4",
        "CaCl2",
        "CaCO3",
        "MgCl2",
        "ZnBr2",
        "DMF",
        "DMSO",
        "THF",
        "DCM",
        "CHCl3",
        "CCl4",
        "MeCN",
        "Acetone",
        "EtOAc",
        "Hexane",
        "Toluene",
        "Benzene",
        "Pyridine",
        "Imidazole",
        "DBU",
        "DABCO",
        "HMPA",
        "HMDS",
        "LiHMDS",
        "NaHMDS",
        "KHMDS",
        "PPh3",
        "P(OEt)3",
        "DEAD",
        "DIAD",
        "AIBN",
        "TEMPO",
        "NaOCl",
        "NaOBr",
        "NaI",
        "KI",
        "CsF",
        "CsCO3",
        "Cs2CO3",
        "AgNO3",
        "Ag2O",
        "AgOAc",
        "Cu(OAc)2",
        "CuSO4",
        "FeSO4",
        "Fe(acac)3",
        "Zn(OAc)2",
        "ZnO",
        "ZnSO4",
        "MnO2",
        "Mn(OAc)2",
        "NiCl2",
        "Ni(acac)2",
        "CoCl2",
        "Co(OAc)2",
        "RuCl3",
        "Ru(bpy)3",
        "IrCl3",
        "PtCl2",
        "PtO2",
        "AuCl",
        "AuCl3",
    ]

    def is_likely_reagent(smiles):
        """Determine if a molecule is likely a reagent rather than a building block."""
        # Check if it's in our common reagents list
        if smiles in common_reagents:
            return True

        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False

            # Very small molecules are likely reagents
            if mol.GetNumAtoms() <= 3:
                return True

            # Check for common reagent functional groups
            if (
                checker.check_fg("Triflate", smiles)
                or checker.check_fg("Mesylate", smiles)
                or checker.check_fg("Tosylate", smiles)
                or checker.check_fg("Acyl halide", smiles)
                or checker.check_fg("Magnesium halide", smiles)
                or checker.check_fg("Zinc halide", smiles)
            ):
                return True

            # Molecules with only C, H, O, N, S, P and <= 8 atoms are often reagents
            atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
            if (
                set(atom_symbols).issubset({"C", "H", "O", "N", "S", "P"})
                and mol.GetNumAtoms() <= 8
            ):
                return True

            # Simple salts are likely reagents
            if (
                any(
                    symbol in ["Li", "Na", "K", "Mg", "Ca", "Cl", "Br", "I", "F"]
                    for symbol in atom_symbols
                )
                and mol.GetNumAtoms() <= 10
            ):
                return True

            return False
        except:
            # If we can't parse it, assume it's not a reagent
            return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early return if we already know it's not linear
            return

        if node["type"] == "reaction":
            # Extract reactants
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Count non-reagent reactants
                building_blocks = []
                for reactant in reactants_smiles:
                    if not is_likely_reagent(reactant):
                        building_blocks.append(reactant)

                if len(building_blocks) > 1:
                    print(
                        f"Detected convergent step at depth {depth} with {len(building_blocks)} major reactants"
                    )
                    is_linear = False
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return is_linear
