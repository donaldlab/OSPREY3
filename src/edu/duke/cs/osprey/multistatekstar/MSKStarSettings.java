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

import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.KStarScore.PartitionFunctionType;
/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
public class MSKStarSettings {
	
	public boolean isReportingProgress;
	public double targetEpsilon;
	public int state;
	public int numTopConfsToSave;
	public MSConfigFileParser cfp;
	public KStarScoreType scoreType;
	public MSSearchProblem[] search;
	public boolean isFinal;
	public LMB[] constraints;
	public PartitionFunctionType[] pfTypes;
	public GMECConfEnergyCalculator.Async[] ecalcs;

	public MSKStarSettings() {}
}
