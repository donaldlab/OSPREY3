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

package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ParamSet;

public class ResidueOrderFactory {
	
	public static ResidueOrder getResidueOrder(ParamSet msParams, MSSearchProblem[][] objFcnSearch) {
		String val = msParams.getValue("RESIDUEORDER");
		switch(val.toLowerCase()) {
		case "staticsequential":
			return new ResidueOrderStaticSequential(objFcnSearch);
		case "staticmindom":
		case "staticobjFunchmean":
		case "dynamicobjfunchmean":
		default:
			throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
		}
	}
	
}
