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

package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;
import java.math.MathContext;

import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.ExpFunction;

public class BoltzmannCalculator {
	
	public static double constRT = PoissonBoltzmannEnergy.constRT;

	public final ExpFunction e;

	public BoltzmannCalculator(MathContext mathContext) {
		e = new ExpFunction(mathContext);
	}
	
	public BigDecimal calc(double energy) {
		return e.exp(-energy/constRT);
	}
}
