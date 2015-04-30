/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.handlempi.MPIMaster;
import edu.duke.cs.osprey.handlempi.MPISlaveTask;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class EnergyMatrixCalculator {
    
    ConfSpace searchSpace;
    ArrayList<Residue> shellResidues;
    
    public EnergyMatrixCalculator(ConfSpace s, ArrayList<Residue> sr) {
        searchSpace = s;
        shellResidues = sr;
    }
    //Calculate a pairwise energy matrix based on a pairwise energy function
    
    
   public EnergyMatrix calcPEM(){
       if(EnvironmentVars.useMPI)
           return calcPEMDistributed();
       else
           return calcPEMLocally();
   }
    
    
    public EnergyMatrix calcPEMLocally(){
        //do the energy calculation here
        
        //MAY WANT TO ARRANGE THIS TO DO EPIC AS WELL
        
        EnergyMatrix ans = new EnergyMatrix(searchSpace);
        
        
        for(int res=0; res<searchSpace.numPos; res++){
            
            TermMinECalculator oneBodyECalc = new TermMinECalculator(searchSpace,shellResidues,false,false,res);
            Object oneBodyE = oneBodyECalc.doCalculation();
            ans.oneBody.set( res, (ArrayList<Double>)oneBodyE );

            for(int res2=0; res2<res; res2++){
                TermMinECalculator pairECalc = new TermMinECalculator(searchSpace,shellResidues,false,false,res,res2);
                Object pairE = pairECalc.doCalculation();
                ans.pairwise.get(res).set( res2, (ArrayList<ArrayList<Double>>)pairE );
            }
        }
        
        return ans;
    }
    
    
    public EnergyMatrix calcPEMDistributed(){
        //do energy calculation on slave nodes via MPI
        
        EnergyMatrix ans = new EnergyMatrix(searchSpace);
        
        MPIMaster mm = MPIMaster.getInstance();//we'll only be running one MPI at once
        ArrayList<MPISlaveTask> tasks = new ArrayList<>();
        
        //generate TermMinECalc objects, in the same order as for local calculation,
        //but this time pass them off to MPI
        for(int res=0; res<searchSpace.numPos; res++){
            
            tasks.add( new TermMinECalculator(searchSpace,shellResidues,false,false,res) );

            for(int res2=0; res2<res; res2++)
                tasks.add( new TermMinECalculator(searchSpace,shellResidues,false,false,res,res2) );
        }
        
        ArrayList<Object> calcResults = mm.handleTasks(tasks);
        
        //Now go through our task results in the same order and put the energies in our matrix
        int resultCount = 0;
        
        for(int res=0; res<searchSpace.numPos; res++){
            
            ans.oneBody.set( res, (ArrayList<Double>)calcResults.get(resultCount) );
            resultCount++;

            for(int res2=0; res2<res; res2++){
                ans.pairwise.get(res).set( res2, (ArrayList<ArrayList<Double>>)calcResults.get(resultCount) );
                resultCount++;
            }
        }
        
        return ans;
    }
    
    
}
