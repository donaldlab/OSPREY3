/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author mhall44
 */
public class ResTemplateMatching {
    
    //the residue and template we're matching
    ResidueTemplate template;
    Residue res;
    
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
    
    public ResTemplateMatching(Residue res, ResidueTemplate template){
        //compute the best matching
        
        this.res = res;
        this.template = template;
        
        score = Double.POSITIVE_INFINITY;//start without any good matching found
        
        ArrayList<Atom> templateAtoms = template.templateRes.atoms;
        numAtoms = templateAtoms.size();
        
        if( res.atoms.size() != numAtoms ){//matching impossible: not even number of atoms matches
            //System.out.println("Missing atoms... Matching will be best, not exact.");
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
        else//we actually only need estimates of bond lengths in the template...other distances don't matter
            //estimate these as the equilibrium bond lengths
            templateDistanceMatrix = template.templateRes.estBondDistanceMatrix();
        
        residueDistanceMatrix = res.atomDistanceMatrix();
        
        
        //if this default-order depth-first search is too slow, maybe try searching
        //least common atoms first (to reduce number of options) or something?
        searchForMatching(0);
    }
    
    
    void searchForMatching(int level){
        //DFS over permutations of atoms, starting at current level (indicates which atom
        //in the template we want to match next)
        
        if(level==numAtoms || level >= res.atoms.size()){//all atoms defined...let's score this matching
            //when we get a full permutation, we score it based on least-square length difference for the template bonds.
            double curMatchingScore = scoreCurPartialMatching();
            if(curMatchingScore<score){//best so far
                score = curMatchingScore;
                System.arraycopy(partialMatching,0,matching,0,numAtoms);//store the matching as well as the score
            }
        }
        else {//need to search more levels

            for(int resAtNum=0; resAtNum<res.atoms.size(); resAtNum++){
                if(!atomAlreadyInPartialMatching(resAtNum,level)){
                    
                    //make sure element types match
                    String templateEleType = template.templateRes.atoms.get(level).elementType;
                    String resEleType = res.atoms.get(resAtNum).elementType;
                    
                    if(templateEleType.equalsIgnoreCase(resEleType)){
                    
                        partialMatching[level] = resAtNum;
                        if(!canPrunePartialMatching(level))//we prune branches if atom pairs that are supposed to be bonded are >2x the distance in the templates
                            searchForMatching(level+1);
                    }
                }
            }
        }
        
    }
    
    boolean atomAlreadyInPartialMatching(int resAtNum, int level){
        //indicates that atom can't be put into partialMatching (currently filling in given level)
        //because it's assigned already at a lower level
        for(int prevLevel=0; prevLevel<level; prevLevel++){
            if(partialMatching[prevLevel] == resAtNum)
                return true;
        }
        return false;//not found
    }
    
    
    boolean canPrunePartialMatching(int level){
        //We just assigned an atom to the given level of the partialMatching
        //see if can prune now.  Assuming can't prune based on prior levels,
        //so will just check bond lengths for newly assigned atom
        
        int newAssignedAtom = partialMatching[level];
        //index of the just-assigned atom in the non-template residue
        
        //level is the index of the just-assigned atom in the template
        for( int bondedTemplateAtNum : templateBonds.get(level) ){
            if(bondedTemplateAtNum<level){//and thus bondedTemplateAtom is already assigned in matching
                //so make sure the bond length is short enough
                
                int bondedAssignedAtNum = partialMatching[bondedTemplateAtNum];
                double templateDist = templateDistanceMatrix[level][bondedTemplateAtNum];
                double resDist = residueDistanceMatrix[newAssignedAtom][bondedAssignedAtNum];
                
                if(resDist >= 1.5*templateDist)//over 1.5 x template length--not a feasible bond
                    //FORMERLY 2 BUT TOO EXPENSIVE
                    return true;
            }
        }
        
        //if we get here, we couldn't prune
        return false;
    }
    
    
    double scoreCurPartialMatching(){
        //score the partial matching based on matches in inter-atom distances
        //for residues that are bonded in the template
        
        double ans = 0;
        
        for(int templateAtNum=0; templateAtNum<numAtoms; templateAtNum++){
            for( int bondedTemplateAtNum : templateBonds.get(templateAtNum) ){
                //we need to compare the templateAtNum-bondedTemplateAtom distance
                //to the distance between the corresponding non-template atoms
                
                int resAtNum = partialMatching[templateAtNum];
                int resBondedAtNum = partialMatching[bondedTemplateAtNum];
                
                double templateDist = templateDistanceMatrix[templateAtNum][bondedTemplateAtNum];
                double resDist = 0;

                if(resAtNum >= 0 && resBondedAtNum >= 0)
                    resDist = residueDistanceMatrix[resAtNum][resBondedAtNum];
                
                double distDiff = templateDist - resDist;
                ans += distDiff*distDiff;
            }
        }
        
        return ans;
        
        //This way we don't have to make assumptions about max bond length (Except for liberal 2x)
        //These assumptions sometimes caused trouble in OSPREY 2
    }
    
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
    
    
}
