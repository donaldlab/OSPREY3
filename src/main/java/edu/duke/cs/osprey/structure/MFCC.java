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

package edu.duke.cs.osprey.structure;


import edu.duke.cs.osprey.tools.PeriodicTable;

import java.util.ArrayList;
import java.util.*;

/**
 *  This class fragments a molecule into residual fragments with conjugate caps
 * @author Hunter Stephens
 */
public class MFCC  {

    public Molecule mol;                     //parent molecule
    public ArrayList<Molecule> fragments;    //capped fragmets
    public ArrayList<Molecule> concaps;      //conjugate caps

    public MFCC(Molecule mol) {

        // set molecule
        this.mol = mol;

        // fractionate and cap
        this.fractionate();
    }

    public void fractionate(){
        int Nresidues = this.mol.residues.size();

        for( int i=1; i < Nresidues-1; i++){
            // create new molecule for fragment
            Molecule frag = new Molecule();

            // get adjacent residues
            frag.appendResidue(this.mol.residues.get(i - 1));
            frag.appendResidue(this.mol.residues.get(i));
            frag.appendResidue(this.mol.residues.get(i+1));

            // take care of dangling bonds
            // Case 1: Left cap is amino terminus
            if(i == 1)
            {
                // replace right terminal C=O  with H

                // 1) find CB and O in fragment
                int Oid = -1, CBid = -1;
                for(int j=0; j < frag.residues.get(2).atoms.size(); j++)
                {
                    if(frag.residues.get(2).atoms.get(j).name.equals("O"))
                    {
                        Oid = j;
                    };
                    if(frag.residues.get(2).atoms.get(j).name.equals("CB"))
                    {
                        CBid = j;
                    };
                }

                // 2) replace CB with H
                frag.residues.get(2).atoms.get(CBid).name = "HleftCap";
                PeriodicTable.setElementProperties(frag.residues.get(0).atoms.get(CBid),"H");

                // 2) cut O
                frag.residues.get(2).atoms.remove((Oid));

            }
            else if(i > 1 && i < Nresidues-1)
            {
                // replace left terminal NH with H
                // replace right terminal C=O with H
            }
            else{
                //replace left terminal NH with H
            }


        }
    }
}
