/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/


package edu.duke.cs.osprey.gmec;

import java.io.Serializable;

/**
 *
 * @author mhall44
 */
public class PruningSettings implements Serializable {
    public boolean typedep = false;
    public double boundsThresh = 100;
    public int algOption = 1;
    public boolean useFlags = true;
    public boolean useTriples = false;
    public double stericThresh = 100;
}

