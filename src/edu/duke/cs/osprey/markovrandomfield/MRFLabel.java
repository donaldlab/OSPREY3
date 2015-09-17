/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.markovrandomfield;

import edu.duke.cs.osprey.confspace.SuperRC;
/**
 *
 * @author hmn5
 */
public class MRFLabel {
    //A label, is an assignment to a node. In our case this will be a super-RC
    //Since the ematrix only cares about the super-RC number, that is how we will store it
    
    //Identifier of label, given node
    int labelNum;
    //the current belief (probability) of the label;
    double currentBelief;
    
    public MRFLabel(int labelNum){
        this.labelNum = labelNum;
    }

}
