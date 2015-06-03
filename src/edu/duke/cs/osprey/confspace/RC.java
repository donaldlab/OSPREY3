/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import java.io.Serializable;
import edu.duke.cs.osprey.tools.MinVolEllipse;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 * @author am439
 */
public class RC implements Serializable {
    //A residue conformation.  Meant to be part of a PositionConfSpace
    
    public String AAType;//amino-acid type
    public int rotNum;//rotamer number (same library, for the given AAType)
    
    //bounds on degrees of freedom
    //some of these are defined by the AAType and rotNum
    public ArrayList<DegreeOfFreedom> DOFs;
    public ArrayList<Double> DOFmin;//minimum values of the DOFs for conformations in the RC
    public ArrayList<Double> DOFmax;
    //note: for AA type we do not use DOFmin or DOFmax (can leave at 0 or whatever): use AAType instead
    
    
    int RCIndex;//index within the RCs for this residue in the PositionConfSpace

    public RC(String AAType, int rotNum, ArrayList<DegreeOfFreedom> DOFs, ArrayList<Double> DOFmin, ArrayList<Double> DOFmax, int RCIndex) {
        this.AAType = AAType;
        this.rotNum = rotNum;
        this.DOFs = DOFs;
        this.DOFmin = DOFmin;
        this.DOFmax = DOFmax;
        this.RCIndex = RCIndex;
    }    
}
