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
