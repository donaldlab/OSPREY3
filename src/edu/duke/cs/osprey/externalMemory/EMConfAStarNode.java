package edu.duke.cs.osprey.externalMemory;

import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;

public class EMConfAStarNode implements ConfAStarNode {
	
	public static final int NotAssigned = -1;
	
	// NOTE: we'll only keep a few nodes in memory at once,
	// so there's less pressure here to optimize space.
	// also, nodes have to be completely self-contained
	// meaning, they can't refer to other nodes for information (eg, links)
	private int[] assignments;
	private int level;
	private double gscore;
	private double hscore;
	
	public EMConfAStarNode(int numPos) {
		assignments = new int[numPos];
		Arrays.fill(assignments, NotAssigned);
        printAssignments();
		level = 0;
		gscore = 0.0;
		hscore = 0.0;
	}
	
	public EMConfAStarNode(EMConfAStarNode other) {
		assignments = Arrays.copyOf(other.assignments, other.assignments.length);
        checkAssignments();
        printAssignments();
		level = other.level;
		gscore = other.gscore;
		hscore = other.hscore;
	}

    public void checkNode()
    {
        System.out.println("Verifying node:");
        printAssignments();
        checkAssignments();
    }

    public String toReadableString()
    {
        String out = this+":\n";
        out+= "level "+level+"\n";
        out+="Assignment:(";
        for(int i = 0; i < assignments.length; i++)
        {
            out+=i+":"+assignments[i];
            if(i < assignments.length-1)
                out+=" ";
        }
        out+=")";
        return out;
    }

    private void printAssignments()
    {
        System.out.println(toReadableString());
    }

    private void checkAssignments()
    {
        for(int i : assignments)
        {
            if(i < -1)
            {
                System.out.println("Error! Offending node:");
                printAssignments();
                throw new Error("negative assignment...");
            }
        }
    }

	@Override
	public EMConfAStarNode assign(int pos, int rc) {
        System.out.println("Assigning "+rc+" to  pos "+pos+"!");
        checkAssignments();
		EMConfAStarNode other = new EMConfAStarNode(this);
        other.checkAssignments();
		other.assignments[pos] = rc;
        if(rc < -1)
            throw new Error("negative assignment...");
		other.level++;
        printAssignments();
        System.out.println("Other: ");
        other.printAssignments();
		return other;
	}

	@Override
	public double getGScore() {
		return gscore;
	}

	@Override
	public void setGScore(double val) {
		gscore = val;
	}

	@Override
	public double getHScore() {
		return hscore;
	}

	@Override
	public void setHScore(double val) {
		hscore = val;
	}

	@Override
	public int getLevel() {
		return level;
	}
	
	public void setLevel(int val) {
        checkAssignments();
        System.out.println("Setting level to "+val);

		level = val;
	}

	@Override
	public void getConf(int[] out) {
		System.arraycopy(assignments, 0, out, 0, assignments.length);
	}
	
	public int[] getConf() {
		return assignments;
	}
	
	@Override
	public void index(ConfIndex index) {
		
        checkAssignments();
		index.node = this;
		
		// copy values and references to the stack for speed
		int n = index.numPos;
		int numDefined = 0;
		int[] dpos = index.definedPos;
		int[] rcs = index.definedRCs;
		int numUndefined = 0;
		int[] upos = index.undefinedPos;
		
		// split conformation into defined and undefined positions
		numDefined = 0;
		numUndefined = 0;
		for (int pos=0; pos<n; pos++) {
			int rc = assignments[pos];
			
			if (rc == EMConfAStarNode.NotAssigned) {
				upos[numUndefined] = pos;
				numUndefined++;
			} else {
				dpos[numDefined] = pos;
                if(rc < 0)
                    throw new Error("Negative RC assignment to ConfIndex position "
                        +pos+":"+rc);
				rcs[numDefined] = rc;
				numDefined++;
			}
		}
		
		// copy stack vars back to the index
		index.numDefined = numDefined;
		index.numUndefined = numUndefined;
	}
}
