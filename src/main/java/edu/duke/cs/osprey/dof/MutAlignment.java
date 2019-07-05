/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.dof;


/* Represents an alignment between two residue templates for use in mutation
With this class we should be able to handle any residue change with the following conditions:
1. There's a "CA-equivalent" atom A
2. There's a "sidechain" whose atoms are only connected to other residues through A
3. The two templates have identical element types and bonds except for the sidechain ("backbone" atoms)
We'll assume here for convenience that these corresponding atoms have the same names
4. There are at least two additional backbone atoms bonded to A, not all collinear
Then the mutation consists of replacing the sidechain atoms with those of the new template,
aligned based on the frame of reference of the three backbone atoms, while the backbone atoms
do not move.

Note mutating to/from proline does not strictly fit this pattern due to the ring;
it's supported as a special case

Also note: if all atoms are equivalent and not standard AAs, this MutAlignment will call everything backbone
and hence mutations will do nothing. (We can't even tell what the sidechain's supposed to be in that case)
*/

import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.InterResBondingTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.tools.VectorAlgebra;

import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class MutAlignment {

    ResidueTemplate oldTemplate, newTemplate;

    TreeSet<String> nonmovingAtomNames;//nice to have a consistent sorted order.  This is the "backbone"
    int[][] mutAlignmentAtoms;//Atom to use to compute the rigid-body motion that puts the "sidechain" in place

    public MutAlignment(ResidueTemplate oldTemplate, ResidueTemplate newTemplate){
        this.oldTemplate = oldTemplate;
        this.newTemplate = newTemplate;
        String curCA = oldTemplate.CAEquivalent;

        if(curCA==null || newTemplate.CAEquivalent==null){
            throw new RuntimeException("ERROR: Trying to align templates " + oldTemplate.name + " and " + newTemplate.name +
                    " for mutation but CAEQUIVALENT not specified in template and cannot be inferred");
        }

        if(!newTemplate.CAEquivalent.equalsIgnoreCase(curCA)) {
            throw makeError("sidechains are defined differently: CA equivalent is "+
                curCA+" for old template, "+newTemplate.CAEquivalent+" for new");
        }
        if(curCA==null)
            throw makeError("sidechain undefined (CAEquivalent==null)");
        if(newTemplate.templateRes.coords==null)//old template doesn't necessarily need coordinates
            throw makeError("new template has no coordinates");

        checkInterResBondingMatch();

        //Pro is a special case: sidechain loops back to N as well as CA
        //Cyx similar because SG eventually bonds back to the backbone
        if(oldTemplate.name.equalsIgnoreCase("PRO")
                || newTemplate.name.equalsIgnoreCase("PRO")
                || oldTemplate.name.equalsIgnoreCase("CYX")
                || newTemplate.name.equalsIgnoreCase("CYX")){
            initStdAAAlignment();//Will detect Pro, handle if possible
            return;
        }
        else {
            matchAtoms();
            handlePeptideTermini();
            checkCAValidity(curCA);

            //kick out atoms from matchedAtoms if they're beyond the CA
            ArrayList<String> matchedSCAtoms = new ArrayList<>();
            for (String at : nonmovingAtomNames) {
                if (getOldTemplateAtom(at) != null && getNewTemplateAtom(at) != null) {
                    if (isMoreSidechainward(at, curCA, true) && isMoreSidechainward(at, curCA, false))
                        matchedSCAtoms.add(at);
                }
            }
            nonmovingAtomNames.removeAll(matchedSCAtoms);

            //at this point can check if the standard N, CA, C preference applies
            if (curCA.equalsIgnoreCase("CA") && nonmovingAtomNames.contains("N") && nonmovingAtomNames.contains("C")) {
                initStdAAAlignment();
                return;
            }

            //get the rest of the reference frame (equivalents to N and C)
            String N = pickRefFrameAtom(curCA, null);
            String C = pickRefFrameAtom(curCA, N);

            String refFrameAtoms[] = new String[]{curCA, N, C};

            mutAlignmentAtoms = new int[2][3];
            for (int a = 0; a < 3; a++) {
                mutAlignmentAtoms[0][a] = oldTemplate.templateRes.getAtomIndexByName(refFrameAtoms[a]);
                mutAlignmentAtoms[1][a] = newTemplate.templateRes.getAtomIndexByName(refFrameAtoms[a]);
            }
        }
    }



    private boolean newTemplateHasAtom(String name){
        //Given the name of an atom in the old template,
        //check that the new template has an equivalent atom (in terms of naming, bonds to the
        //current atoms in nonMovingAtoms, and element)
        Atom oldTemplateAtom = getOldTemplateAtom(name);
        Atom newTemplateAtom = getNewTemplateAtom(name);
        if(newTemplateAtom==null)
            return false;
        if(newTemplateAtom.elementNumber != oldTemplateAtom.elementNumber)
            return false;
        //now check bonding equivalence
        Predicate<String> nm = a->nonmovingAtomNames.contains(a);
        Set<String> oldNonMovingBonds = oldTemplateAtom.bonds.stream().map(b->b.name).filter(nm).collect(Collectors.toSet());
        Set<String> newNonMovingBonds = newTemplateAtom.bonds.stream().map(b->b.name).filter(nm).collect(Collectors.toSet());
        for(String bondedName : oldNonMovingBonds){
            if(!newNonMovingBonds.contains(bondedName))
                return false;
        }
        for(String bondedName : newNonMovingBonds){
            if(!oldNonMovingBonds.contains(bondedName))
                return false;
        }
        return true;
    }

    private void matchAtoms(){
        //Match atoms between the old and new templates, and store their names in nonmovingAtomNames
        nonmovingAtomNames = new TreeSet<>();
        for (Atom at : oldTemplate.templateRes.atoms) {
            if (oldTemplate.interResBonding.atomCanBondOtherRes(at)) {
                if (newTemplateHasAtom(at.name))
                    nonmovingAtomNames.add(at.name);
                else
                    throw makeError("atom " + at.name + " bonded to other residue but changes in mutation");
            }
        }

        if (nonmovingAtomNames.isEmpty()) {
            //if no inter-residue bonds, just start from the first atom that appears in both old & new template
            for (Atom at : oldTemplate.templateRes.atoms) {
                if (newTemplateHasAtom(at.name)) {
                    nonmovingAtomNames.add(at.name);
                    break;
                }
            }
        }

        if (nonmovingAtomNames.isEmpty()) {
            throw new RuntimeException("ERROR: Can't mutate " + oldTemplate.name + " to "
                    + newTemplate.name + " because they have no atoms in common");
        }

        //OK now expand the list of matched atoms outward until no longer equivalent
        while(true) {
            Atom nextAtom = null;
            for (String prevMatchedAtom : nonmovingAtomNames) {
                for (Atom at2 : getOldTemplateAtom(prevMatchedAtom).bonds) {
                    if (newTemplateHasAtom(at2.name) && !nonmovingAtomNames.contains(at2.name)) {
                        nextAtom = at2;
                        break;
                    }
                }
                if (nextAtom != null)
                    break;
            }
            if(nextAtom==null)
                break;
            else
                nonmovingAtomNames.add(nextAtom.name);
        }
    }


    private boolean canBeCA(String possibleCA){
        //Can the atom with the specified name be used as the root of a mutating sidechain?
        //i.e. as a CA equivalent
        for(Atom at : oldTemplate.templateRes.atoms){
            if(!nonmovingAtomNames.contains(at.name)) {//at is supposed to be a sidechain atom
                if (!isMoreSidechainward(at.name, possibleCA, true)){
                    return false;//can't be CA because at is not on a branch off of possibleCA
                }
            }
        }
        for(Atom at : newTemplate.templateRes.atoms){
            if(!nonmovingAtomNames.contains(at.name)) {//at is supposed to be a sidechain atom
                if (!isMoreSidechainward(at.name, possibleCA, false)){
                    return false;//can't be CA because at is not on a branch off of possibleCA
                }
            }
        }
        return true;
    }

    private boolean isMoreSidechainward(String at1, String at2, boolean useOldTemplate){
        //at1 is more sidechainward than at2 if at1 can only reach other residues via at2
        HashSet<String> atomsVisited = new HashSet<>();
        atomsVisited.add(at2);
        return ! canEscapeResNotVia(at1, atomsVisited, useOldTemplate);
    }

    private boolean canEscapeResNotVia(String atName, HashSet<String> forbidden, boolean useOldTemplate){
        //Is there a path of bonds from at to another res that doesn't go through the forbidden atoms?
        ResidueTemplate templ = useOldTemplate ? oldTemplate : newTemplate;
        Atom at = templ.templateRes.getAtomByName(atName);
        if(templ.interResBonding.atomCanBondOtherRes(at))
            return true;
        forbidden.add(atName);//don't loop back indefinitely, in case of cycles
        for(Atom at2 : at.bonds){
            if(!forbidden.contains(at2.name)) {
                if (canEscapeResNotVia(at2.name, forbidden, useOldTemplate))
                    return true;
            }
        }
        return false;
    }


    private String pickRefFrameAtom(String CA, String N){
        //We've picked a CA, now pick another reference frame atom
        //If we already have one ("N"!=null), get another that's not collinear
        //prefer atoms bonded to CA (to get a more rigid reference frame)
        //note the "CA" and "N" equivalents don't actually have to be carbon and nitrogen
        //I'm just calling them that by analogy to (alpha) amino acids
        for(Atom cand : getOldTemplateAtom(CA).bonds){
            if(nonmovingAtomNames.contains(cand.name)){
                if(isGoodRefFrameAtom(cand.name,CA,N))
                    return cand.name;
            }
        }
        for(String candName : nonmovingAtomNames){
            if(isGoodRefFrameAtom(candName,CA,N))
                return candName;
        }

        throw makeError("can't find 3 noncollinear unchanging (backbone) atoms to use as a frame of reference");
    }


    private boolean isGoodRefFrameAtom(String atom, String CA, String N) {
        //Is this atom good for defining the reference frame for mutation?
        //The CA-equivalent and, if not null, the N-equivalent are provided
        if(!nonmovingAtomNames.contains(atom))
            return false;
        if(getNewTemplateAtom(atom)==null || getOldTemplateAtom(atom)==null)
            return false;
        if(atom.equalsIgnoreCase(N) || atom.equalsIgnoreCase(CA))//already taken
            return false;
        if(N!=null) {//make sure the atoms are non-collinear
            //and in approximately the same geometry between the two templates
            //we can check these things by looking at the distances between the atoms
            ArrayList<Double> newDists = interAtomDistances(new String[] {atom,CA,N}, false);
            if(oldTemplate.templateRes.coords!=null) {//we may not always have coords for the old template (need them for new one though)
                ArrayList<Double> oldDists = interAtomDistances(new String[]{atom, CA, N}, true);
                for (int a = 0; a < 3; a++) {
                    if (Math.abs(oldDists.get(a) - newDists.get(a)) > 0.1 * Math.abs(oldDists.get(a)))
                        return false;//check the triangles are congruent (or very close)
                }
            }
            if(Collections.max(newDists)>0.499*(newDists.get(0)+newDists.get(1)+newDists.get(2)))
                return false;//corresponds to at least ~170 degree angle
        }
        return true;
    }

    private Atom getOldTemplateAtom(String name) {
        return oldTemplate.templateRes.getAtomByName(name);
    }

    private Atom getNewTemplateAtom(String name) {
        return newTemplate.templateRes.getAtomByName(name);
    }

    private void initStdAAAlignment(){
        if( ! (HardCodedResidueInfo.hasAminoAcidBB(oldTemplate.templateRes)&&HardCodedResidueInfo.hasAminoAcidBB(newTemplate.templateRes)) ){
            throw makeError(" can't do standard AA alignment for non-amino acids");
        }
        nonmovingAtomNames = new TreeSet<>(HardCodedResidueInfo.listBBAtomsForMut(newTemplate,oldTemplate));
        mutAlignmentAtoms = HardCodedResidueInfo.findMutAlignmentAtoms(oldTemplate,newTemplate);
    }

    private void handlePeptideTermini(){
        //special case: allow mutating away charged peptide termini
        //this means putting the backbone H's and O's in the backbone even if they disappear on mutation
        if( HardCodedResidueInfo.hasAminoAcidBB(oldTemplate.templateRes)
                && HardCodedResidueInfo.hasAminoAcidBB(newTemplate.templateRes) ) {
            Atom N = getOldTemplateAtom("N");
            Atom C = getOldTemplateAtom("C");
            if (oldTemplate.interResBonding.atomCanBondOtherRes(N)
                    && oldTemplate.interResBonding.atomCanBondOtherRes(C)) {
                //no sidechain atoms bond N or C except in Pro, which has already been pulled aside
                //by the time we call this function
                //Also deal with newTemplate: if old has H1, H2, H3 new may have H, etc.
                for(Atom nc : new Atom[] {N,C,getNewTemplateAtom("N"),getNewTemplateAtom("C")}){
                    for (Atom at : nc.bonds)
                        nonmovingAtomNames.add(at.name);
                }
            }
        }
    }

    private static Set<String> getBondingAtomNames(InterResBondingTemplate irb){
        if(irb instanceof InterResBondingTemplate.SpecifiedBondingAtomsTemplate)
            return new HashSet<>(((InterResBondingTemplate.SpecifiedBondingAtomsTemplate) irb).getBondingAtomNames());
        else if(irb instanceof InterResBondingTemplate.NoBondingTemplate)
            return new HashSet<>();
        else if(irb instanceof InterResBondingTemplate.CysteineBondingTemplate)
            return new HashSet<>(Arrays.asList("N","C","SG"));
        else if(irb instanceof InterResBondingTemplate.PeptideBondingTemplate)
            return new HashSet<>(Arrays.asList("N","C"));
        else
            throw new RuntimeException("ERROR: Unknown inter-residue bonding type");
    }

    private void checkInterResBondingMatch(){
        //complain if residue have different inter-res bonding, preventing mutation
        InterResBondingTemplate irb1 = oldTemplate.interResBonding;
        InterResBondingTemplate irb2 = newTemplate.interResBonding;
        if(irb1 instanceof InterResBondingTemplate.PeptideBondingTemplate
                && irb2 instanceof InterResBondingTemplate.PeptideBondingTemplate)
            return;//OK (even if one is cysteine)
        else {
            Set<String> bondingAtoms1 = getBondingAtomNames(irb1);
            Set<String> bondingAtoms2 = getBondingAtomNames(irb2);
            for(String name : bondingAtoms1){
                if(!bondingAtoms2.contains(name))
                    throw makeError("Only the old template has inter-residue bonding through atom "+name);
            }
            for(String name : bondingAtoms2){
                if(!bondingAtoms1.contains(name))
                    throw makeError("Only the new template has inter-residue bonding through atom "+name);
            }
            return;
        }
    }

    private ArrayList<Double> interAtomDistances(String[] atomNames, boolean useOldTemplate){
        Function<String,Atom> at = s -> (useOldTemplate ? getOldTemplateAtom(s) : getNewTemplateAtom(s));
        List<Atom> atoms = Arrays.asList(atomNames).stream().map(at).collect(Collectors.toList());
        ArrayList<Double> ans = new ArrayList<>();
        ans.add(VectorAlgebra.distance(atoms.get(0).getCoords(), atoms.get(1).getCoords()));
        ans.add(VectorAlgebra.distance(atoms.get(2).getCoords(), atoms.get(1).getCoords()));
        ans.add(VectorAlgebra.distance(atoms.get(0).getCoords(), atoms.get(2).getCoords()));
        return ans;
    }

    private boolean areThereMovingAtoms(){
        for(ResidueTemplate templ : new ResidueTemplate[] {oldTemplate,newTemplate}){
            for(Atom at : templ.templateRes.atoms){
                if(!nonmovingAtomNames.contains(at.name))
                    return true;
            }
        }
        return false;
    }

    private String inferCA(){
        //check that the non-nonMovingAtoms form a valid "sidechain", and find its CA
        //(a single atom that they all go through)
        //NOT using this now: inferring on a per-mutation basis lets us use an arbitrary sidechain
        //for mutating original structure to template for original res type; thus no coordinates change
        //and this leads to inconsistency because eventually we'll mutate to something else and back
        //and the template could be subtly but noticeably different from the original crystal structure
        boolean CAisCA = false;//there is an atom called CA that can function as the CA
        String curCA = null;
        for (String at : nonmovingAtomNames) {
            if(getNewTemplateAtom(at)!=null && getOldTemplateAtom(at)!=null) {//CA must be truly shared, not special peptide terminal atom
                if (canBeCA(at)) {
                    if (at.equalsIgnoreCase("CA"))
                        CAisCA = true;
                    if (curCA == null)
                        curCA = at;
                    else if (isMoreSidechainward(curCA, at, true) && isMoreSidechainward(curCA, at, false))
                        curCA = at;
                    //let's choose the CA option closest to the other res, as in standard AA's. In case geometry around CB-equivalent changes
                }
            }
        }

        if(curCA==null)
            throw makeError("changing atoms are not all in a single sidechain");

        return curCA;
    }

    private void checkCAValidity(String curCA){
        if(getOldTemplateAtom(curCA)==null || getNewTemplateAtom(curCA)==null)
            throw makeError("invalid CA equivalent: "+curCA+".  Must appear in both old and new template");
        if(!canBeCA(curCA))
            throw makeError("invalid CA equivalent: "+curCA+".  Some changing atoms are not in the sidechain starting at it");
    }


    private RuntimeException makeError(String specifMessage){
        //Make an error indicating we can't do this mutation, and why
        return new RuntimeException("ERROR: Can't mutate "+oldTemplate.name+" to "+newTemplate.name+" because "+specifMessage);
    }

    public int[][] getMutAlignmentAtoms(){
        return mutAlignmentAtoms;
    }

    public TreeSet<String> getNonmovingAtomNames(){
        return nonmovingAtomNames;
    }
}
