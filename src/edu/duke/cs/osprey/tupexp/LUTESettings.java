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


package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.control.ParamSet;
import java.io.Serializable;

/**
 *
 * Settings for LUTE
 * 
 * @author mhall44
 */
@SuppressWarnings("serial")
public class LUTESettings implements Serializable {
    
    boolean useLUTE=false;
    public double goalResid=0.01;
    boolean useRelWt=false;
    boolean useThreshWt=false;
    
    public LUTESettings(){
        //by default, no LUTE
        useLUTE = false;
    }
    
    public LUTESettings(ParamSet params){
        //initialize from input parameter set
        useLUTE = params.getBool("USETUPEXP");
        goalResid = params.getDouble("LUTEGOALRESID");
        useRelWt = params.getBool("LUTERELWT");
        useRelWt = params.getBool("LUTETHRESHWT");
    }
    
    
    public static LUTESettings defaultLUTE(){
        LUTESettings ans = new LUTESettings();
        ans.useLUTE = true;
        return ans;
    }
    
    public boolean shouldWeUseLUTE(){
        return useLUTE;
    }
}

