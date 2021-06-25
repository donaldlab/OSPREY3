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

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import edu.duke.cs.osprey.control.Main;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import org.apache.commons.lang3.text.WordUtils;

import edu.duke.cs.osprey.structure.Residue.SecondaryStructure;
import edu.duke.cs.osprey.tools.FileTools;


/**
 * XYZ format class for Quantum Chem input files
 *
 * @author Hunter Stephens
 */
public class XYZIO {

    public static String write(Molecule mol){
        StringBuilder buf = new StringBuilder();

        for (Residue res : mol.residues) {
            for (Atom atom : res.atoms) {
                buf.append(atom.elementType);
                buf.append("\t");
                for(int i=0; i<3; i++){
                    buf.append(atom.getCoords()[i]);
                    buf.append("\t");
                }
                buf.append("\n");
            }
        }
        return buf.toString();
    }
}
