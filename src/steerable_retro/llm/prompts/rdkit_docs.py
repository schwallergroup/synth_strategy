docs = """class rdkit.Chem.rdFMCS.AtomCompare
CompareAny = rdkit.Chem.rdFMCS.AtomCompare.CompareAny
CompareAnyHeavyAtom = rdkit.Chem.rdFMCS.AtomCompare.CompareAnyHeavyAtom
CompareElements = rdkit.Chem.rdFMCS.AtomCompare.CompareElements
CompareIsotopes = rdkit.Chem.rdFMCS.AtomCompare.CompareIsotopes

class rdkit.Chem.rdFMCS.BondCompare
CompareAny = rdkit.Chem.rdFMCS.BondCompare.CompareAny
CompareOrder = rdkit.Chem.rdFMCS.BondCompare.CompareOrder
CompareOrderExact = rdkit.Chem.rdFMCS.BondCompare.CompareOrderExact

rdkit.Chem.rdFMCS.FindMCS((AtomPairsParameters)mols[, (bool)maximizeBonds=True[, (float)threshold=1.0[, (int)timeout=3600[, (bool)verbose=False[, (bool)matchValences=False[, (bool)ringMatchesRingOnly=False[, (bool)completeRingsOnly=False[, (bool)matchChiralTag=False[, (AtomCompare)atomCompare=rdkit.Chem.rdFMCS.AtomCompare.CompareElements[, (BondCompare)bondCompare=rdkit.Chem.rdFMCS.BondCompare.CompareOrder[, (RingCompare)ringCompare=rdkit.Chem.rdFMCS.RingCompare.IgnoreRingFusion[, (str)seedSmarts='']]]]]]]]]]]]) → MCSResult
Find the MCS for a set of molecules

class rdkit.Chem.rdFMCS.MCSAcceptance((object)arg1)
Base class. Subclass and override MCSAcceptance.__call__() to define a custom boolean callback function. Returning True will cause the MCS candidate to be accepted, False to be rejected

class rdkit.Chem.rdFMCS.MCSAtomCompare((object)arg1)
Base class. Subclass and override MCSAtomCompare.__call__() to define custom atom compare functions, then set MCSParameters.AtomTyper to an instance of the subclass

CheckAtomCharge((MCSAtomCompare)self, (MCSAtomCompareParameters)parameters, (Mol)mol1, (int)atom1, (Mol)mol2, (int)atom2) → bool
Return True if both atoms have the same formal charge

CheckAtomChirality((MCSAtomCompare)self, (MCSAtomCompareParameters)parameters, (Mol)mol1, (int)atom1, (Mol)mol2, (int)atom2) → bool
Return True if both atoms have, or have not, a chiral tag

CheckAtomRingMatch((MCSAtomCompare)self, (MCSAtomCompareParameters)parameters, (Mol)mol1, (int)atom1, (Mol)mol2, (int)atom2) → bool
Return True if both atoms are, or are not, in a ring

class rdkit.Chem.rdFMCS.MCSAtomCompareParameters((object)arg1)
Parameters controlling how atom-atom matching is done

class rdkit.Chem.rdFMCS.MCSBondCompare((object)arg1)
Base class. Subclass and override MCSBondCompare.__call__() to define custom bond compare functions, then set MCSParameters.BondTyper to an instance of the subclass

CheckBondRingMatch((MCSBondCompare)self, (MCSBondCompareParameters)parameters, (Mol)mol1, (int)bond1, (Mol)mol2, (int)bond2) → bool
Return True if both bonds are, or are not, part of a ring

CheckBondStereo((MCSBondCompare)self, (MCSBondCompareParameters)parameters, (Mol)mol1, (int)bond1, (Mol)mol2, (int)bond2) → bool
Return True if both bonds have, or have not, a stereo descriptor

class rdkit.Chem.rdFMCS.MCSBondCompareParameters((object)arg1)
Parameters controlling how bond-bond matching is done

class rdkit.Chem.rdFMCS.MCSFinalMatchCheck((object)arg1)
Base class. Subclass and override MCSFinalMatchCheck.__call__() to define a custom boolean callback function. Returning True will cause the growing seed to be accepted, False to be rejected

class rdkit.Chem.rdFMCS.MCSParameters((object)arg1)
Parameters controlling how the MCS is constructed

class rdkit.Chem.rdFMCS.MCSProgress((object)arg1)
Base class. Subclass and override MCSProgress.__call__() to define a custom callback function

class rdkit.Chem.rdFMCS.MCSProgressData((object)arg1)
Information about the MCS progress

class rdkit.Chem.rdFMCS.MCSResult
Bases: instance

used to return MCS results

Raises an exception This class cannot be instantiated from Python

property canceled
if True, the MCS calculation did not finish
property degenerateSmartsQueryMolDict
Dictionary collecting all degenerate (SMARTS, queryMol) pairs (empty if MCSParameters.StoreAll is False)
property numAtoms
number of atoms in MCS
property numBonds
number of bonds in MCS
property queryMol
query molecule for the MCS
property smartsString
SMARTS string for the MCS

class rdkit.Chem.rdFMCS.RingCompare
###
IgnoreRingFusion = rdkit.Chem.rdFMCS.RingCompare.IgnoreRingFusion
PermissiveRingFusion = rdkit.Chem.rdFMCS.RingCompare.PermissiveRingFusion
StrictRingFusion = rdkit.Chem.rdFMCS.RingCompare.StrictRingFusion
###
class rdkit.Chem.rdchem.RingInfo
RingInfo objects store information about the rings in a molecule.
###

The following functions are available for the RingInfo class:
    AddRing((RingInfo)self, (AtomPairsParameters)atomIds, (AtomPairsParameters)bondIds) → None
    Adds a ring to the set. Be very careful with this operation.

    AreAtomsInSameRing((RingInfo)self, (int)idx1, (int)idx2) → bool

    AreAtomsInSameRingOfSize((RingInfo)self, (int)idx1, (int)idx2, (int)size) → bool

    AreBondsInSameRing((RingInfo)self, (int)idx1, (int)idx2) → bool

    AreBondsInSameRingOfSize((RingInfo)self, (int)idx1, (int)idx2, (int)size) → bool

    AreRingFamiliesInitialized((RingInfo)self) → bool

    AreRingsFused((RingInfo)self, (int)ring1Idx, (int)ring2Idx) → bool

    AtomMembers((RingInfo)self, (int)idx) → object

    AtomRingFamilies((RingInfo)self) → object

    AtomRingSizes((RingInfo)self, (int)idx) → object

    AtomRings((RingInfo)self) → object

    BondMembers((RingInfo)self, (int)idx) → object

    BondRingFamilies((RingInfo)self) → object

    BondRingSizes((RingInfo)self, (int)idx) → object

    BondRings((RingInfo)self) → object

    IsAtomInRingOfSize((RingInfo)self, (int)idx, (int)size) → bool

    IsBondInRingOfSize((RingInfo)self, (int)idx, (int)size) → bool

    IsRingFused((RingInfo)self, (int)ringIdx) → bool

    MinAtomRingSize((RingInfo)self, (int)idx) → int

    MinBondRingSize((RingInfo)self, (int)idx) → int

    NumAtomRings((RingInfo)self, (int)idx) → int

    NumBondRings((RingInfo)self, (int)idx) → int

    NumFusedBonds((RingInfo)self, (int)ringIdx) → int

    NumRelevantCycles((RingInfo)self) → int

    NumRingFamilies((RingInfo)self) → int

    NumRings((RingInfo)self) → int
###
class rdkit.Chem.rdchem.Mol((object)self)
###

The following functions are available for the Mol class:

    AddConformer((Mol)self, (Conformer)conf[, (bool)assignId=False]) → int
    Add a conformer to the molecule and return the conformer ID

    ClearComputedProps((Mol)self[, (bool)includeRings=True]) → None
    Removes all computed properties from the molecule.

    ClearProp((Mol)self, (str)key) → None
    Removes a property from the molecule.

    Compute2DCoords((Mol)mol[, (bool)canonOrient=True[, (bool)clearConfs=True[, (dict)coordMap={}[, (int)nFlipsPerSample=0[, (int)nSample=0[, (int)sampleSeed=0[, (bool)permuteDeg4Nodes=False[, (float)bondLength=-1.0[, (bool)forceRDKit=False[, (bool)useRingTemplates=False]]]]]]]]]]) → int
    Compute 2D coordinates for a molecule.

    ComputeGasteigerCharges((Mol)mol[, (int)nIter=12[, (bool)throwOnParamFailure=False]]) → None
    Compute Gasteiger partial charges for molecule

    Debug(useStdout=False)

    GetAromaticAtoms((Mol)self) → _ROQAtomSeq
    Returns a read-only sequence containing all of the molecule’s aromatic Atoms.

    GetAtomWithIdx((Mol)self, (int)idx) → Atom
    Returns a particular Atom.

    GetAtoms()

    GetAtomsMatchingQuery((Mol)self, (QueryAtom)qa) → _ROQAtomSeq
    Returns a read-only sequence containing all of the atoms in a molecule that match the query atom.

    GetBondBetweenAtoms((Mol)self, (int)idx1, (int)idx2) → Bond
    Returns the bond between two atoms, if there is one.

    GetBondWithIdx((Mol)self, (int)idx) → Bond
    Returns a particular Bond.

    GetBonds()

    GetBoolProp((Mol)self, (str)key) → bool
    Returns the Bool value of the property if possible.

    GetConformer((Mol)self[, (int)id=-1]) → Conformer
    Get the conformer with a specified ID

    GetConformers((Mol)self) → _ROConformerSeq
    Returns a read-only sequence containing all of the molecule’s Conformers.

    GetDoubleProp((Mol)self, (str)key) → float
    Returns the double value of the property if possible.

    GetIntProp((Mol)self, (str)key) → int
    Returns the integer value of the property if possible.

    GetNumAtoms((Mol)self[, (int)onlyHeavy=-1[, (bool)onlyExplicit=True]]) → int
    Returns the number of atoms in the molecule.

    GetNumBonds((Mol)self[, (bool)onlyHeavy=True]) → int
    Returns the number of Bonds in the molecule.

    GetNumConformers((Mol)self) → int
    Return the number of conformations on the molecule

    GetNumHeavyAtoms((Mol)self) → int
    Returns the number of heavy atoms (atomic number >1) in the molecule.

    GetProp((Mol)self, (str)key[, (bool)autoConvert=False]) → object
    Returns the value of the property.

    GetPropNames((Mol)self[, (bool)includePrivate=False[, (bool)includeComputed=False]]) → _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
    Returns a tuple with all property names for this molecule.

    GetPropsAsDict((Mol)self[, (bool)includePrivate=False[, (bool)includeComputed=False[, (bool)autoConvertStrings=True]]]) → dict
    Returns a dictionary populated with the molecules properties.

    GetRingInfo((Mol)self) → RingInfo
    Returns the number of molecule’s RingInfo object.

    GetStereoGroups((Mol)self) → StereoGroup_vect
    Returns a list of StereoGroups defining the relative stereochemistry of the atoms.

    GetSubstructMatch((Mol)self, (Mol)query[, (bool)useChirality=False[, (bool)useQueryQueryMatches=False]]) → object
    Returns the indices of the molecule’s atoms that match a substructure query.

    GetSubstructMatches((Mol)self, (Mol)query[, (bool)uniquify=True[, (bool)useChirality=False[, (bool)useQueryQueryMatches=False[, (int)maxMatches=1000]]]]) → object
    Returns tuples of the indices of the molecule’s atoms that match a substructure query.

    GetUnsignedProp((Mol)self, (str)key) → int
    Returns the unsigned int value of the property if possible.

    HasProp((Mol)self, (str)key) → int
    Queries a molecule to see if a particular property has been assigned.

    HasQuery((Mol)self) → bool
    Returns if any atom or bond in molecule has a query

    HasSubstructMatch((Mol)self, (Mol)query[, (bool)recursionPossible=True[, (bool)useChirality=False[, (bool)useQueryQueryMatches=False]]]) → bool
    Queries whether or not the molecule contains a particular substructure.

    NeedsUpdatePropertyCache((Mol)self) → bool
    Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.

    RemoveAllConformers((Mol)self) → None
    Remove all the conformations on the molecule

    RemoveConformer((Mol)self, (int)id) → None
    Remove the conformer with the specified ID

    SetBoolProp((Mol)self, (str)key, (bool)val[, (bool)computed=False]) → None
    Sets a boolean valued molecular property

    SetDoubleProp((Mol)self, (str)key, (float)val[, (bool)computed=False]) → None
    Sets a double valued molecular property

    SetIntProp((Mol)self, (str)key, (int)val[, (bool)computed=False]) → None
    Sets an integer valued molecular property

    SetProp((Mol)self, (str)key, (str)val[, (bool)computed=False]) → None
    Sets a molecular property

    SetUnsignedProp((Mol)self, (str)key, (int)val[, (bool)computed=False]) → None
    Sets an unsigned integer valued molecular property

    ToBinary((Mol)self) → object
    Returns a binary string representation of the molecule.

    UpdatePropertyCache((Mol)self[, (bool)strict=True]) → None
    Regenerates computed properties like implicit valence and ring information.
###
class rdkit.Chem.rdChemReactions.ChemicalReaction((object)self)
###

The following functions are available for the ChemicalReaction class:
    AddAgentTemplate((ChemicalReaction)self, (Mol)mol) → int
    adds a agent (a Molecule)

    AddProductTemplate((ChemicalReaction)self, (Mol)mol) → int
    adds a product (a Molecule)

    AddReactantTemplate((ChemicalReaction)self, (Mol)mol) → int
    adds a reactant (a Molecule) to the reaction

    AddRecursiveQueriesToReaction((ChemicalReaction)reaction[, (dict)queries={}[, (str)propName='molFileValue'[, (bool)getLabels=False]]]) → object
    adds recursive queries and returns reactant labels

    ClearComputedProps((ChemicalReaction)self) → None
    Removes all computed properties from the reaction.

    ClearProp((ChemicalReaction)self, (str)key) → None
    Removes a property from the reaction.

    GetAgentTemplate((ChemicalReaction)self, (int)which) → Mol
    returns one of our agent templates

    GetAgents((ChemicalReaction)self) → MOL_SPTR_VECT
    get the agent templates

    GetBoolProp((ChemicalReaction)self, (str)key) → bool
    Returns the Bool value of the property if possible.

    GetDoubleProp((ChemicalReaction)self, (str)key) → float
    Returns the double value of the property if possible.

    GetIntProp((ChemicalReaction)self, (str)key) → int
    Returns the integer value of the property if possible.

    GetNumAgentTemplates((ChemicalReaction)self) → int
    returns the number of agents this reaction expects

    GetNumProductTemplates((ChemicalReaction)self) → int
    returns the number of products this reaction generates

    GetNumReactantTemplates((ChemicalReaction)self) → int
    returns the number of reactants this reaction expects

    GetProductTemplate((ChemicalReaction)self, (int)which) → Mol
    returns one of our product templates

    GetProducts((ChemicalReaction)self) → MOL_SPTR_VECT
    get the product templates

    GetProp((ChemicalReaction)self, (str)key) → str
    Returns the value of the property.

    GetPropNames((ChemicalReaction)self[, (bool)includePrivate=False[, (bool)includeComputed=False]]) → _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
    Returns a tuple with all property names for this reaction.

    GetPropsAsDict((ChemicalReaction)self[, (bool)includePrivate=False[, (bool)includeComputed=False[, (bool)autoConvertStrings=True]]]) → dict
    Returns a dictionary populated with the reaction’s properties.

    GetReactantTemplate((ChemicalReaction)self, (int)which) → Mol
    returns one of our reactant templates

    GetReactants((ChemicalReaction)self) → MOL_SPTR_VECT
    get the reactant templates

    GetReactingAtoms((ChemicalReaction)self[, (bool)mappedAtomsOnly=False]) → object
    returns a sequence of sequences with the atoms that change in the reaction

    GetSubstructParams((ChemicalReaction)self) → SubstructMatchParameters
    get the parameter object controlling the substructure matching

    GetUnsignedProp((ChemicalReaction)self, (str)key) → int
    Returns the unsigned int value of the property if possible.

    HasProp((ChemicalReaction)self, (str)key) → int
    Queries a molecule to see if a particular property has been assigned.

    Initialize((ChemicalReaction)self[, (bool)silent=False]) → None
    initializes the reaction so that it can be used

    IsInitialized((ChemicalReaction)self) → bool
    checks if the reaction is ready for use

    IsMoleculeAgent((ChemicalReaction)self, (Mol)mol) → bool
    returns whether or not the molecule has a substructure match to one of the agents.

    IsMoleculeProduct((ChemicalReaction)self, (Mol)mol) → bool
    returns whether or not the molecule has a substructure match to one of the products.

    IsMoleculeReactant((ChemicalReaction)self, (Mol)mol) → bool
    returns whether or not the molecule has a substructure match to one of the reactants.

    RemoveAgentTemplates((ChemicalReaction)self[, (AtomPairsParameters)targetList=None]) → None
    Removes agents from reaction. If targetList is provide the agents will be transferred to that list.

    RemoveUnmappedProductTemplates((ChemicalReaction)self[, (float)thresholdUnmappedAtoms=0.2[, (bool)moveToAgentTemplates=True[, (AtomPairsParameters)targetList=None]]]) → None
    Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from product templates to the agent templates or to a given targetList

    RemoveUnmappedReactantTemplates((ChemicalReaction)self[, (float)thresholdUnmappedAtoms=0.2[, (bool)moveToAgentTemplates=True[, (AtomPairsParameters)targetList=None]]]) → None
    Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from reactant templates to the agent templates or to a given targetList

    RunReactant((ChemicalReaction)self, (AtomPairsParameters)reactant, (int)reactionIdx) → object
    apply the reaction to a single reactant

    RunReactantInPlace((ChemicalReaction)self, (Mol)reactant[, (bool)removeUnmatchedAtoms=True]) → bool
    apply the reaction to a single reactant in place. The reactant itself is modified. This can only be used for single reactant - single product reactions.

    RunReactants((ChemicalReaction)self, (tuple)reactants[, (int)maxProducts=1000]) → object
    apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples. If maxProducts is not zero, stop the reaction when maxProducts have been generated [default=1000]

    SetBoolProp((ChemicalReaction)self, (str)key, (bool)val[, (bool)computed=False]) → None
    Sets a boolean valued molecular property

    SetDoubleProp((ChemicalReaction)self, (str)key, (float)val[, (bool)computed=False]) → None
    Sets a double valued molecular property

    SetIntProp((ChemicalReaction)self, (str)key, (int)val[, (bool)computed=False]) → None
    Sets an integer valued molecular property

    SetProp((ChemicalReaction)self, (str)key, (str)val[, (bool)computed=False]) → None
    Sets a molecular property

    SetUnsignedProp((ChemicalReaction)self, (str)key, (int)val[, (bool)computed=False]) → None
    Sets an unsigned integer valued molecular property

    ToBinary((ChemicalReaction)self) → object
    Returns a binary string representation of the reaction.

    Validate((ChemicalReaction)self[, (bool)silent=False]) → tuple
    checks the reaction for potential problems, returns (numWarnings,numErrors)
###
"""