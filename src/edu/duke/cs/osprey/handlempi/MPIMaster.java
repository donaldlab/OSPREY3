/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.handlempi;

import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class MPIMaster {
    //This is an object to handle MPI master functions
    //The program has just one of these, so it's a singleton class OR JUST DO EVERYTHING STATIC
    //When given a list of MPISlaveTasks to handle, it'll farm them out to slave nodes
    //and return their calculation results
    
    
    static int processRank = 0;
    
    private MPIMaster() {
    }
    
    public static void printIfMaster(String output){
        //print the string if called at the master node.  Otherwise do nothing.  
        if(processRank==0)
            System.out.println(output);
    }
    
    public static MPIMaster getInstance() {
        if(processRank>0)
            throw new RuntimeException("ERROR: Slave nodes can't get the MPIMaster instance");
        
        return MPIMasterHolder.INSTANCE;
    }

    public ArrayList<Object> handleTasks(ArrayList<MPISlaveTask> tasks) {
        //Given the list of tasks, return their results in the same order
        
        //DEBUG!!!!
        //BEFORE WE IMPLEMENT MPI JUST DO THESE LOCALLY
        ArrayList<Object> ans = new ArrayList<>();
        for(MPISlaveTask task : tasks){
            ans.add( task.doCalculation() );
        }
        return ans;
    }
    
    private static class MPIMasterHolder {

        private static final MPIMaster INSTANCE = new MPIMaster();
    }
}
