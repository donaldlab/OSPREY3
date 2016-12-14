/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper;

/**
 * 
 * This class checks backbone dihedrals against experimental Ramachandran-plot contours.
 * 
 * @author mhall44
 */

import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import java.io.*;
import java.util.Arrays;
import java.util.StringTokenizer;



public class RamachandranChecker {


    double[][][] tables;//tables[a][b][c] is density at phi=b+/-1, psi=c+/-1
    //with a=0 is for gly, a=1 for pro, a=2 for general, and a=3 for pre-pro

//This class is designed to read the Richardsons' top500 Ramachandran plot density data:
//(from top500-angles/pct/rama)

    //2-degree bins for each angle, centered at -179, -177,..
    //Density (0-to-1 scale) in 3rd column; sorted by phi and then by psi;
    //columns separated by spaces; starts with some comment lines marked with a #



    static double denCutoff = 0.02f;//Cutoff density for being allowed

    private RamachandranChecker() {
        
    }

    public static RamachandranChecker getInstance() {
        return RamachandranCheckerHolder.INSTANCE;
    }

    private static class RamachandranCheckerHolder {
        private static final RamachandranChecker INSTANCE = new RamachandranChecker();
    }


    public void readInputFiles(String[] fileNames){//Filenames in order: gly, pro, general, pre-pro

        tables = new double[4][180][180];

        for(int a=0;a<4;a++){

            try{
                BufferedReader br = new BufferedReader( new FileReader( fileNames[a] ) );
                String line = br.readLine();

                while( line.charAt(0) == '#' )
                    line = br.readLine();

                for(int phiBin=0; phiBin<180; phiBin++){
                    for(int psiBin=0; psiBin<180; psiBin++){

                        StringTokenizer st = new StringTokenizer(line," ");

                        st.nextToken();
                        st.nextToken();

                        tables[a][phiBin][psiBin] = Double.valueOf(st.nextToken());

                        line = br.readLine();
                    }
                }

                br.close();
            }
            catch(IOException ex){
                throw new Error("can't read Ramachandran plot file " + fileNames[a], ex);
            }

        }
    }
    
    
    public boolean[] checkByAAType(Residue res){
        //Return the acceptability at {gly, pro, other AA types} of a given residue's BB dihedrals
        
        boolean ans[] = new boolean[3];

        double phiPsi[] = getPhiPsi(res);
        
        if(phiPsi==null){//undefined
            Arrays.fill(ans, true);//can't rule conformation out!
            return ans;
        }

        for(int a=0;a<3;a++)
            ans[a] = checkAngles(phiPsi[0], phiPsi[1], a);

        return ans;
    }


    //Same but for prePro
    public boolean checkPrePro(Residue res){

        double phiPsi[] = getPhiPsi(res);

        if(phiPsi==null)//undefined
            return true;
                    
        return checkAngles(phiPsi[0], phiPsi[1], 3);
    }

    
    //Returns {phi,psi} for the residue. 
    public double[] getPhiPsi(Residue res){

        double ans[] = new double[2];

        //Get coordinates of relevant atoms
        //return null (undefined) if can't find one or more atoms
        if(res.indexInMolecule==0 || res.indexInMolecule==res.molec.residues.size()-1)//first or last res
            return null;
        
        Residue prevRes = res.molec.residues.get(res.indexInMolecule-1);
        Residue nextRes = res.molec.residues.get(res.indexInMolecule+1);
        
        
        double[] CLast = prevRes.getCoordsByAtomName("C");
        double[] NCur = res.getCoordsByAtomName("N");
        double[] CACur = res.getCoordsByAtomName("CA");
        double[] CCur = res.getCoordsByAtomName("C");
        double[] NNext = nextRes.getCoordsByAtomName("N");
        
        if ( CLast==null || NCur==null || CACur==null || CCur==null || NNext==null )
            return null;//atom not found
        
        ans[0] = Protractor.measureDihedral( new double[][] {CLast,NCur,CACur,CCur} );//phi
        ans[1] = Protractor.measureDihedral( new double[][] {NCur,CACur,CCur,NNext} );//psi

        return ans;
    }



    public boolean checkAngles(double phi, double psi, int plotNum){
        
        phi = getInRange(phi);
        psi = getInRange(psi);
        
        int phiBin = (int)((phi+180)/2);
        int psiBin = (int)((psi+180)/2);
        double den = tables[plotNum][phiBin][psiBin];
        if(den > denCutoff)
            return true;
        else
            return false;
    }

    
    double getInRange(double angle){
        //get angle in the range [-180,180), which is required for both phi and psi
        while(angle>=180)
            angle -= 360;
        while(angle<-180)
            angle += 360;
        return angle;
    }


 }

