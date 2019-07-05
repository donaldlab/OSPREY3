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

package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;

import java.util.*;
import java.util.stream.Collectors;

/**
 *
 * @author mhall44
 */
public class ResTemplateMatching {
    
    //the residue and template we're matching
    ResidueTemplate template;
    Residue res;
    ForcefieldParams ffParams;
    
    int[] matching;//for each atom in the template, index of the corresponding atom in the residue
    
    int[] partialMatching;//We'll use this during our depth-first search for matchings;
    //any good ones will be copied to matching
    
    public double score;//least-square deviations in bond lengths, for atoms pairs bonded in the template
    //if all atoms bonded in the template are bonded, then we expect this is the right template and matching
    
    int numAtoms;//number of atoms in the residue (should be the same as in the template)
    
    
    //some things useful in performing the assignment
    ArrayList<ArrayList<Integer>> templateBonds;//For each template atom, list of indices of other atoms
    //in the template that it's bonded to
    double[][] residueDistanceMatrix;//matrix of distances between all atoms in the non-template residue
    //(ordered as in the residue)
    double[][] templateDistanceMatrix;
    
    
    //for dynamically ordered DFBB
    int[] templateAtomOrdering;//At each level we assign template atom templateAtomOrdering[level]
    int[] templateAtomOrderingRev;//reverse mapping of templateOrdering
    double cumulativeScores[];//Cache of the score contributions from each level of our search
    
    boolean use13Distances = false;//Use (1,3) bond distances for matching as well
    //as direct bonds.  This had some issues with proline (C and CB seem weirdly close?)
    //So will have to be used with caution
    
    int curAssignmentCount = 0;//how many atom assignments we've tried so far
    
    static double timeoutableScore = 0.01;//if taking a long time, this is the per-
    //atom score that is considered good enough to stop
    int timeoutNumAssignments = 100000000;
    
    public ResTemplateMatching(Residue res, ResidueTemplate template, ForcefieldParams ffParams){
        //compute the best matching
        
        this.res = res;
        this.template = template;
        this.ffParams = ffParams;
        
        score = Double.POSITIVE_INFINITY;//start without any good matching found
        
        ArrayList<Atom> templateAtoms = template.templateRes.atoms;
        numAtoms = templateAtoms.size();
        
        if( res.atoms.size() != numAtoms ){//matching impossible: not even number of atoms matches
            return;
        }
        
        //initialize search
        matching = new int[templateAtoms.size()];
        
        partialMatching = new int[templateAtoms.size()];
        Arrays.fill(partialMatching,-1);//all undefined now
        
        //compute the helper quantities
        templateBonds = new ArrayList<>();
        for(Atom atom1 : template.templateRes.atoms){
            ArrayList<Integer> atom1Bonds = new ArrayList<>();
            for(int atNum2=0; atNum2<numAtoms; atNum2++){
                Atom atom2 = template.templateRes.atoms.get(atNum2);
                for(Atom bondedAtom : atom1.bonds){
                    if(bondedAtom==atom2){
                        atom1Bonds.add(atNum2);
                        break;
                    }
                }
            }
            templateBonds.add(atom1Bonds);
        }
        
        if(template.templateRes.coords!=null)//we have coordinates for templateRes...can use real coordinates
            templateDistanceMatrix = template.templateRes.atomDistanceMatrix();
        else {//we actually only need estimates of bond lengths in the template...other distances don't matter
            //estimate these as the equilibrium bond lengths
            templateDistanceMatrix = template.templateRes.estBondDistanceMatrix(ffParams);
            use13Distances = false;
        }
        
        residueDistanceMatrix = res.atomDistanceMatrix();
        
        chooseTemplateAtomOrdering();
        cumulativeScores = new double[numAtoms];//so far no scores yet
        
        //if this default-order depth-first search is too slow, maybe try searching
        //least common atoms first (to reduce number of options) or something?
        //Also can prune much better by using partial scores
        //can have an array storing the cumulative score up to each level
        //levels only add to score so it's a good lower bound for pruning
        //can test against full score at the end
        //can test with 1AK0.pdb.  bew_bes.bin is there already so comment template making if desired
        curAssignmentCount = 0;
        searchForMatching(0);
    }
    
    
    private void searchForMatching(int level){
        //DFBB over permutations of atoms, starting at current level
        //next template atom we want to match is templateAtomOrdering[level]
        //any time we find the best matching so far, store it in matching
        
        if(curAssignmentCount>timeoutNumAssignments){
            //taking a while...see if a good enough score has been found
            if(score < timeoutableScore * res.atoms.size())
                return;
        }
        
        double curMatchingScore = scoreCurPartialMatching(level);
        
        if(level==numAtoms){//all atoms defined
            if(curMatchingScore<score){//best so far
                if(checkStereochemistry()) {//check we didn't mess up the stereochemistry before accepting the matching
                    score = curMatchingScore;
                    System.arraycopy(partialMatching, 0, matching, 0, numAtoms);//store the matching as well as the score
                }
            }
        }
        else if(curMatchingScore<score){
            //need to search more levels unless current score is already no better than
            //our current best score (score can only rise as we add levels)  
            //note this criterion will remove any impossible (infinite-scored) matchings
            
            //If the residue has an atom matching the template atom's name, try it first
            String templateAtomName = template.templateRes.atoms.
                    get(templateAtomOrdering[level]).name;
            int nameMatchAtom = res.getAtomIndexByName(templateAtomName);
            if(nameMatchAtom>-1)
                tryNextAtomMatch(level, nameMatchAtom);
            
            for(int resAtNum=0; resAtNum<numAtoms; resAtNum++){//number of atom in the residue we're trying to match
                if(resAtNum!=nameMatchAtom)
                    tryNextAtomMatch(level, resAtNum);
            }
        }
    }
    
    private void tryNextAtomMatch(int level, int resAtNum){
        //see if template atom templateAtomOrdering[level] can match residue atom resAtNum
        curAssignmentCount++;
        if(!atomAlreadyInPartialMatching(resAtNum,level)){
            //make sure element types match
            int templateEleType = template.templateRes.atoms.
                    get(templateAtomOrdering[level]).elementNumber;
            int resEleType = res.atoms.get(resAtNum).elementNumber;

            if(templateEleType == resEleType){
                partialMatching[templateAtomOrdering[level]] = resAtNum;
                searchForMatching(level+1);
            }
        }
    }
    
    boolean atomAlreadyInPartialMatching(int resAtNum, int level){
        //indicates that atom can't be put into partialMatching (currently filling in given level)
        //because it's assigned already at a lower level
        for(int prevLevel=0; prevLevel<level; prevLevel++){
            if(partialMatching[templateAtomOrdering[prevLevel]] == resAtNum)
                return true;
        }
        return false;//not found
    }
    
    Collection<Integer> getConstrDistAtoms(int atomNum1){
        //Given the index in template of an atom 1, get the indices in template
        //of other atoms that are expected to be at a fixed distance from atom 1
        if(use13Distances){
            HashSet<Integer> ans = new HashSet<>();
            for(int bondedAtNum : templateBonds.get(atomNum1)){
                ans.add(bondedAtNum);
                for(int atNum13 : templateBonds.get(bondedAtNum)){
                    if(atNum13!=atomNum1)
                        ans.add(atNum13);
                }
            }
            return ans;
        }
        else
            return templateBonds.get(atomNum1);
        //note: by only considering distances between atoms that we are assigning
        //as bonded or 1,3-bonded, we avoid artificially creating bonds when there are clashes
        //between atoms that the template says are not bonded
    }
    
    double scoreCurPartialMatching(int level){
        //Score the current partial matching: sum over all pairs of bonded
        //and assigned template atoms of the bond distance deviation
        //all deviations not involving the current atom are cached in cumulativeScores
        
        //first score the most recent assignment (level-1).  Cache the score terms
        if(level>0){
            cumulativeScores[level-1] = 0;
            int justAssignedTemplateAtom = templateAtomOrdering[level-1];
            int justAssignedResAtom = partialMatching[justAssignedTemplateAtom];
            for( int bondedTemplateAtNum : getConstrDistAtoms(justAssignedTemplateAtom) ){
                if(templateAtomOrderingRev[bondedTemplateAtNum]<level){//bondedTemplateAtom already assigned in partialMatching
                    
                    int bondedResAtNum = partialMatching[bondedTemplateAtNum];
                    double templateDist = templateDistanceMatrix[justAssignedTemplateAtom][bondedTemplateAtNum];
                    double resDist = residueDistanceMatrix[justAssignedResAtom][bondedResAtNum];

                    if(resDist >= 1.5*templateDist){//over 1.5 x template length--not a feasible bond
                        cumulativeScores[level-1] = Double.POSITIVE_INFINITY;//prune this horrible conformation
                        break;
                    }
                    else {
                        double distDiff = templateDist - resDist;
                        cumulativeScores[level-1] += distDiff*distDiff;
                    }
                }
            }
        }
            
        //Every assignedLevel<level can contribute to the score
        double ans = 0;
        for(int assignedLevel=0; assignedLevel<level; assignedLevel++){
            ans += cumulativeScores[assignedLevel];
        }
        return ans;
    }
    
    //This way we don't have to make assumptions about max bond length (Except for liberal 1.5x)
    //unless we need to generate a new template, in which case we can have more user supervision
    //These assumptions sometimes caused trouble in OSPREY 2
    
    
    public void assign(){
        //assign the template

        res.template = template;
        
        //copy atoms from template
        ArrayList<Atom> newAtoms = new ArrayList<>();
        for(int atNum=0; atNum<numAtoms; atNum++){
            Atom newAtom = template.templateRes.atoms.get(atNum).copy();
            newAtom.res = res;
            newAtoms.add(newAtom);
        }
        res.atoms = newAtoms;
        res.markIntraResBondsByTemplate();
        
        //reorder coordinates to match the template ordering of atoms
        double[] newCoords = new double[3*numAtoms];
        for(int atNum=0; atNum<numAtoms; atNum++){
            System.arraycopy(res.coords, 3*matching[atNum], newCoords, 3*atNum, 3);
        }
        res.coords = newCoords;
    }
    
    
    private void chooseTemplateAtomOrdering(){
        //decide what order to search the atoms in
        //It is better to try to match rare atoms first
        //because then they can help with the scoring of the more common atoms,
        //thus reducing the amount of branching
        
        HashMap<Integer,Integer> elementFreqs = new HashMap<>();
        //how many atoms there are of each element type
        
        for(Atom at : template.templateRes.atoms){
            if( elementFreqs.containsKey(at.elementNumber) )
                elementFreqs.put(at.elementNumber, elementFreqs.get(at.elementNumber)+1 );
            else
                elementFreqs.put( at.elementNumber, 1 );
        }
        
        //Get all the unique element frequencies, sorted.  HashSet is for uniqueness
        ArrayList<Integer> elementFreqsSorted = 
                new ArrayList<>(new HashSet<>(elementFreqs.values()));
        Collections.sort(elementFreqsSorted);
        
        //Now put the template's atom numbers in array ordered by their element frequency
        int count = 0;
        templateAtomOrdering = new int[numAtoms];
        templateAtomOrderingRev = new int[numAtoms];
        for(int elementFreq : elementFreqsSorted){
            for(int templateAtNum=0; templateAtNum<numAtoms; templateAtNum++){
                int templAtElement = template.templateRes.atoms.get(templateAtNum).elementNumber;
                if(elementFreqs.get(templAtElement)==elementFreq){
                    //put in this atom
                    templateAtomOrdering[count] = templateAtNum;
                    templateAtomOrderingRev[templateAtNum] = count;
                    count++;
                }
            }
        }
        if(count!=numAtoms)
            throw new RuntimeException("ERROR: Bug in ResTemplateMatching, lost "+(numAtoms-count)+" atoms");
    }

    private boolean checkStereochemistry(){
        //Check if the current partial matching flips the stereochemistry at any atom with 4 bonds
        //ONLY WORKS FOR FULLY SPECIFIED PARTIAL MATCHING
        //This fails not only if the template is different stereochemically from the residue,
        //but if any stereochemically distinct groups (e.g. the two CD methyls on leucine) are swapped
        if(template.templateRes.coords==null)//template has no coordinates so we'll just have to assume the stereochemistry is right
            return true;
        for(int atNum=0; atNum<numAtoms; atNum++) {
            Atom templateAtom = template.templateRes.atoms.get(atNum);
            if(templateAtom.bonds.size() == 4){
                Atom resAtom = res.atoms.get(partialMatching[atNum]);
                boolean templateStereoSign = stereoSign(templateAtom, templateAtom.bonds);
                List<Atom> resBondedAtoms = templateAtom.bonds.stream().map(at->res.atoms.get(partialMatching[at.indexInRes])).collect(Collectors.toList());
                boolean resStereoSign = stereoSign(resAtom, resBondedAtoms);
                if(templateStereoSign != resStereoSign)
                    return false;
            }
        }
        return true;
    }

    private static boolean stereoSign(Atom at, List<Atom> subst){
        if(subst.size()!=4)
            throw new RuntimeException("ERROR: don't know how to check stereochemistry with "+subst.size()+" substituents");
        double[] center = at.getCoords();
        List<double[]> bonds = subst.stream().map( bat -> VectorAlgebra.subtract( bat.getCoords(), center ) ).collect(Collectors.toList());
        double[] normal = VectorAlgebra.cross(bonds.get(0), bonds.get(1));
        return VectorAlgebra.dot(normal, bonds.get(2)) > VectorAlgebra.dot(normal, bonds.get(3));
    }

}
