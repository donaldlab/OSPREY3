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

import java.util.ArrayList;

public interface ResidueOrder {

	public class AAScore {
		public int residuePos;
		public int AATypePos;
		public double score;

		public AAScore(int residuePos, int AATypePos, double score) {
			this.residuePos = residuePos;
			this.AATypePos = AATypePos;
			this.score = score;
		}
	}
	
	public ArrayList<ArrayList<ArrayList<AAScore>>> getNextAssignments(MSSearchProblem[][] objFcnSearch, int numMaxMut);
	
}
