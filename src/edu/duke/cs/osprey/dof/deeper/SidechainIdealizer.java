/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper;

import edu.duke.cs.osprey.dof.DihedralRotation;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;

/**
 *
 * This class handles idealization of sidechains;
 * residues moved by DEEPer have idealized sidechains because orig sidechain geometry
 * not valid for perturbed backbones
 * 
 * @author mhall44
 */
public class SidechainIdealizer {
    
    
    
    //Idealize the sidechain of residue res.
    //This function is based on chiropraxis.sc.SidechainIdealizer.idealizeCB from KiNG
    //The whole sidechain idealization (SidechainIdealizer.idealizeSidechain) is not performed because
    //a mutated residue would already be in an idealized conformation, and if the residue is not already
    //(e.g. because it's a wild-type rotamer) it probably shouldn't be
    //This function rotates the sidechain as a rigid body about CA
    public static void idealizeSidechain(ResidueTemplateLibrary templateLib, Residue res){

        //Coordinates of key atoms in the residue
        double NCoord[] = res.getCoordsByAtomName("N");
        double CACoord[] = res.getCoordsByAtomName("CA");
        double CCoord[] = res.getCoordsByAtomName("C");

        if(res.template.name.equalsIgnoreCase("GLY")){
            //No rotation matrix but two HAs to handle

            int HASystem = 0;//The system used to name HA's.  0=HA2/HA3, 1=1HA/2HA, 2=2HA/3HA

            int HANum = res.getAtomIndexByName("HA2");
            if( HANum == -1 ){//It was called 1HA not HA2 maybe
                HANum = res.getAtomIndexByName("1HA");
                HASystem = 1;
            }
            if( HANum == -1 ){//Try the 2HA/3HA system
                HANum = res.getAtomIndexByName("2HA");
                HASystem = 2;
            }
            if( HANum == -1 ){
                throw new Error("could not parse HA names for " + res.fullName);
            }


            double t1[] = VectorAlgebra.get4thPoint(NCoord, CCoord, CACoord, 1.100f, 109.3f, -121.6f);
            double t2[] = VectorAlgebra.get4thPoint(CCoord, NCoord, CACoord, 1.100f, 109.3f, 121.6f);

            double newHA[] = VectorAlgebra.subtract( VectorAlgebra.scale( VectorAlgebra.add(t1,t2) , 0.5f), CACoord);
            newHA = VectorAlgebra.add( VectorAlgebra.scale( newHA, 1.100f / VectorAlgebra.norm(newHA) ), CACoord );

            System.arraycopy(newHA, 0, res.coords, HANum*3, 3);
            
            //now do the CB-like HA
            
            if(HASystem == 0)
                HANum = res.getAtomIndexByName("HA3");
            else if (HASystem == 1)
                HANum = res.getAtomIndexByName("2HA");
            else//HASystem == 2
                HANum = res.getAtomIndexByName("3HA");

            if( HANum == -1 ){
                throw new Error("could not parse HA names for " + res.fullName);
            }

            t1 = VectorAlgebra.get4thPoint(NCoord, CCoord, CACoord, 1.100f, 109.3f, 121.6f);
            t2 = VectorAlgebra.get4thPoint(CCoord, NCoord, CACoord, 1.100f, 109.3f, -121.6f);

            newHA = VectorAlgebra.subtract( VectorAlgebra.scale( VectorAlgebra.add(t1,t2) , 0.5f), CACoord);
            newHA = VectorAlgebra.add( VectorAlgebra.scale( newHA, 1.100f / VectorAlgebra.norm(newHA) ), CACoord );

            System.arraycopy(newHA, 0, res.coords, HANum*3, 3);//Change the HA coordinates
        }
        else{//Will have a C-beta

            double CBCoord[] = res.getCoordsByAtomName("CB");

            double[] t1, t2;

            if(res.template.name.equalsIgnoreCase("ALA")) {
                t1 = VectorAlgebra.get4thPoint(NCoord, CCoord, CACoord, 1.536f, 110.1f, 122.9f);
                t2 = VectorAlgebra.get4thPoint(CCoord, NCoord, CACoord, 1.536f, 110.6f, -122.6f);
            }
            else if(res.template.name.equalsIgnoreCase("PRO")) {
                t1 = VectorAlgebra.get4thPoint(NCoord, CCoord, CACoord, 1.530f, 112.2f, 115.1f);
                t2 = VectorAlgebra.get4thPoint(CCoord, NCoord, CACoord, 1.530f, 103.0f, -120.7f);
            }
            else if( res.template.name.equalsIgnoreCase("VAL") || res.template.name.equalsIgnoreCase("THR") || res.template.name.equalsIgnoreCase("ILE") ) {
                t1 = VectorAlgebra.get4thPoint(NCoord, CCoord, CACoord, 1.540f, 109.1f, 123.4f);
                t2 = VectorAlgebra.get4thPoint(CCoord, NCoord, CACoord, 1.540f, 111.5f, -122.0f);
            }
            else { // anything else
                t1 = VectorAlgebra.get4thPoint(NCoord, CCoord, CACoord, 1.530f, 110.1f, 122.8f);
                t2 = VectorAlgebra.get4thPoint(CCoord, NCoord, CACoord, 1.530f, 110.5f, -122.6f);
            }

            double idealCB[] = VectorAlgebra.scale(VectorAlgebra.add(t1, t2), 0.5f);

            //Now the rotation matrix for the sidechain will give the smallest rotation
            //that maps CBCoord to idealCB
            double ABOld[] = VectorAlgebra.subtract(CBCoord, CACoord);
            double ABNew[] = VectorAlgebra.subtract(idealCB, CACoord);
            double theta = Protractor.getAngleRadians( ABOld, ABNew );
            double rotax[] = VectorAlgebra.cross(ABOld, ABNew);
            //RotationMatrix SC_mtx;

            if( VectorAlgebra.norm(rotax) != 0 ){//CB not already in the ideal position
                RotationMatrix SC_mtx = new RotationMatrix(rotax[0], rotax[1], rotax[2], theta, true);
                RigidBodyMotion SC_motion = new RigidBodyMotion(CACoord, SC_mtx,CACoord);
                moveSidechain(res,SC_motion);
            }

            t1 = VectorAlgebra.get4thPoint(NCoord, CCoord, CACoord, 1.100f, 107.9f, -118.3f);
            t2 = VectorAlgebra.get4thPoint(CCoord, NCoord, CACoord, 1.100f, 108.1f, 118.2f);
            double newHA[] = VectorAlgebra.subtract( VectorAlgebra.scale( VectorAlgebra.add(t1,t2) , 0.5f), CACoord);
            newHA = VectorAlgebra.add( VectorAlgebra.scale( newHA, 1.100f / VectorAlgebra.norm(newHA) ), CACoord );
            
            int HANum = res.getAtomIndexByName("HA");
            System.arraycopy(newHA, 0, res.coords, HANum*3, 3);//Change the HA coordinates
        }

        if(res.template.name.equalsIgnoreCase("PRO"))//Idealize the CG, CD and ring hydrogen positions
            idealizeProRing(templateLib, res);
        //For proline, the above sidechain idealization operations do not move CD
    }
    
    
    
    public static void moveSidechain(Residue res, RigidBodyMotion motion){
        //apply the given motion to the sidechain of res
        int numAtoms = res.atoms.size();
                
        for(int atomNum=0; atomNum<numAtoms; atomNum++){
            String atomName = res.atoms.get(atomNum).name;
            
            if(isSidechainAtom(atomName,res))
                motion.transform(res.coords, atomNum);
        }
    }
    
    
    static boolean isSidechainAtom(String atomName, Residue res){
        //check by name if the atom is in the sidechain (i.e., the atom from CB onward,
        //which are treated together for idealization purposes)
        if(atomName.startsWith("HA"))
            return false;
        if(res.template.name.equalsIgnoreCase("PRO")){
            if(atomName.startsWith("CD")){//Pro CD moved as part of backbone
                return false;
            }
        }

        for(String name2 : HardCodedResidueInfo.possibleBBAtoms){
            if(name2.equalsIgnoreCase(atomName)){
                return false;
            }
        }
        
        return true;
    }
    
    
    
    public static boolean idealizeProRing(ResidueTemplateLibrary templateLib, Residue res){
        //Make an idealized ring for the given residue, given the current positions of the backbone atoms, CB, and CD
        //The pucker specified in res.pucker will be applied:
        //UP pucker means the smaller of the two possible solutions for chi1; DOWN pucker means the larger.
        //The CA-N-CD and CA-CB-CG angles are adjusted to the ideal as they are likely to be too large after mutations, etc.
        //Then chi1 is calculated, finding the two (UP and DOWN) solutions that maintain the ideal CG-CD bond length
        //(from the Amber template)
        //Returns true for success, false for failure

        
        if( ! res.template.name.equalsIgnoreCase("PRO") ){
            throw new Error("trying to idealize the proline ring on a non-proline");
        }
        
        double CANCD_ideal = (double)(Math.PI*111.5/180);
        double CACBCG_ideal = (double)(Math.PI*103.9/180);
        //from Ho et al, Protein Sci 2005;14:1011
        
        Residue idealPro = templateLib.getTemplateForMutation("PRO", res).templateRes;
        

        //actual coordinates
        double NCoord[] = res.getCoordsByAtomName("N");
        double CACoord[] = res.getCoordsByAtomName("CA");
        double CBCoord[] = res.getCoordsByAtomName("CB");
        double CGCoord[] = res.getCoordsByAtomName("CG");
        double CDCoord[] = res.getCoordsByAtomName("CD");


        //First adjust the bond angles
        double ncd[] = VectorAlgebra.subtract( CDCoord, NCoord );
        double canhat[] = VectorAlgebra.subtract( NCoord, CACoord );//Unit vector in CA --> N direction
        double vhat[] = VectorAlgebra.perpendicularComponent(ncd, canhat);
        canhat = VectorAlgebra.scale( canhat, 1/VectorAlgebra.norm(canhat) );
        vhat = VectorAlgebra.scale( vhat, 1/VectorAlgebra.norm(vhat) );
        ncd = VectorAlgebra.scale( VectorAlgebra.add( VectorAlgebra.scale( canhat, -(double)Math.cos(CANCD_ideal) ) ,
                VectorAlgebra.scale( vhat, (double)Math.sin(CANCD_ideal) ) ), VectorAlgebra.norm(ncd) );
        CDCoord = VectorAlgebra.add( ncd , NCoord );

        int CDNum = res.getAtomIndexByName("CD");
        System.arraycopy( CDCoord, 0, res.coords, 3*CDNum, 3 );//Store the new CD coordinates

        double cbcg[] = VectorAlgebra.subtract( CGCoord, CBCoord );
        double cacbhat[] = VectorAlgebra.subtract( CBCoord, CACoord );
        double what[] = VectorAlgebra.perpendicularComponent(cbcg, cacbhat);
        cacbhat = VectorAlgebra.scale( cacbhat, 1/VectorAlgebra.norm(cacbhat) );//Unit vector in CA --> CB direction
        what = VectorAlgebra.scale( what, 1/VectorAlgebra.norm(what) );
        cbcg = VectorAlgebra.scale( VectorAlgebra.add( VectorAlgebra.scale( cacbhat, -(double)Math.cos(CACBCG_ideal) ) ,
                VectorAlgebra.scale( what, (double)Math.sin(CACBCG_ideal) ) ), VectorAlgebra.norm(cbcg) );
        CGCoord = VectorAlgebra.add( cbcg , CBCoord );

        int CGNum = res.getAtomIndexByName("CG");
        System.arraycopy( CGCoord, 0, res.coords, 3*CGNum, 3 );//Store the new CG coordinates (their torsion will be modified later)


        //Now calculate chi1.

        double dir2[] = VectorAlgebra.perpendicularComponent( canhat , cacbhat);
        dir2 = VectorAlgebra.scale( dir2 , 1/VectorAlgebra.norm(dir2) );
        double dir3[] = VectorAlgebra.cross( cacbhat , dir2 );
        double dBG = VectorAlgebra.norm(cbcg);
        double sina = (double)Math.sin(CACBCG_ideal);
        double cosa = (double)Math.cos(CACBCG_ideal);
        //The position of CG will be given by CBCoord + dBG*(-cosa*cacbhat + sina*(dir2*cos(chi1)+dir3*sin(chi1)))

        double CGCD = VectorAlgebra.norm( VectorAlgebra.subtract( idealPro.getCoordsByAtomName("CD") , idealPro.getCoordsByAtomName("CG") ) );
        //We now get chi1 by making the CG-CD distance be chi1

        double[] cbcd = VectorAlgebra.subtract( CDCoord, CBCoord );
        double r1 = VectorAlgebra.dot( cbcd , cacbhat );
        double r2 = VectorAlgebra.dot( cbcd , dir2 );
        double r3 = VectorAlgebra.dot( cbcd , dir3 );

        double A = ( (r1+cosa*dBG)*(r1+cosa*dBG) + r2*r2 + r3*r3 + dBG*dBG*sina*sina - CGCD*CGCD ) / (2*dBG*sina);
        double R = (double)Math.sqrt( r2*r2 + r3*r3 );

        if( Math.abs(A/R) > 1 ){
            //idealizing the ring is impossible
            //mark as a ConfProblem for this res associated with its pucker DOF
            res.pucker.setRingCloseable(false);
            return false;
        }
        else{
            res.pucker.setRingCloseable(true);
            
            double beta = (double)Math.atan2(r2, r3);
            double asr = (double)Math.asin(A/R);

            double newChi1;
            
            //we want to idealize the pucker given by res.pucker.curPucker
            if(res.pucker.getCurPucker() == ProlinePucker.Direction.UP)
                newChi1 = asr - beta;
            else//DOWN
                newChi1 = (double)Math.PI - beta - asr;//Larger newChi1, since an arcsine is acute (puckers merge when asr = pi/2, i.e. A=R)


            newChi1 = 180*newChi1/((double)Math.PI);//Convert to degrees

            //Now move CG into place according to the new chi1
            double dihCoords[][] = new double[][] {NCoord,CACoord,CBCoord,CGCoord};
            double curChi1 = Protractor.measureDihedral(dihCoords);
            DihedralRotation dihRot = new DihedralRotation(CACoord, CBCoord, newChi1-curChi1);
            dihRot.transform(res.coords, CGNum);
            //This will mess up the CG-HG bonds but idealizeRingHydrogens will fix that

            idealizeRingHydrogens(res, idealPro);//Put the hydrogens in place

            return true;
        }
    }


    public static void idealizeRingHydrogens(Residue res, Residue idealPro){
        //Idealize the sidechain hydrogens (those bonded to CB, CG, and CD) of a proline
        //Does not change the heavy-atom coordinates

        double HBAngle = 108.78f*(double)Math.PI/180;//HB1-CB-HB2 angle
        double HGAngle = 108.47f*(double)Math.PI/180;
        double HDAngle = 108.46f*(double)Math.PI/180;
        //These values are from Allen et al, Chem. Eur. J. 2004, 10, 4512 – 4517
        //which is an ab initio study of free neutral proline
        //They provide constraints for structural refinement given heavy-atom positions,
        //including constant values for the above angles plus linear combinations of angles like CA-CB-HB1
        //It is probably not worthwhile to include these small corrections at runtime
        //(seems to require solving a trig equation that's an 8th-order polynomial in half-angle tangent form,
        //plus the polypeptide chain geometry might change things)
        //So we assume C_2v symmetry (i.e. 2-fold rotational symmetry about heavy-atom-angle bisector)
        //at CB, CG, and CD with the above angles

        double HBDist = VectorAlgebra.norm( VectorAlgebra.subtract( idealPro.getCoordsByAtomName("HB2") , idealPro.getCoordsByAtomName("CB") ) );//CB-HB2 (or, assumed, CB-HB1 distance) (HB2 is nice because both HB1/2 and HB2/3 nomenclatures have it)
        double HGDist = VectorAlgebra.norm( VectorAlgebra.subtract( idealPro.getCoordsByAtomName("HG2") , idealPro.getCoordsByAtomName("CG") ) );
        double HDDist = VectorAlgebra.norm( VectorAlgebra.subtract( idealPro.getCoordsByAtomName("HD2") , idealPro.getCoordsByAtomName("CD") ) );


        double CACoord[] = res.getCoordsByAtomName("CA");
        double CBCoord[] = res.getCoordsByAtomName("CB");
        double CGCoord[] = res.getCoordsByAtomName("CG");
        double CDCoord[] = res.getCoordsByAtomName("CD");
        double NCoord[] = res.getCoordsByAtomName("N");

        double CACB[] = VectorAlgebra.subtract(CBCoord,CACoord);
        double CBCG[] = VectorAlgebra.subtract(CGCoord,CBCoord);
        double CGCD[] = VectorAlgebra.subtract(CDCoord,CGCoord);
        double CDN[] = VectorAlgebra.subtract(NCoord,CDCoord);


        double BBisector[] = VectorAlgebra.average( CACB, VectorAlgebra.scale(CBCG,-1) );//Bisector of the HB1-CB-HB2 angle
        double BPerpendicular[] = VectorAlgebra.cross( CACB, CBCG );//Perpendicular to the CA-CB-CG plane
        double GBisector[] = VectorAlgebra.average( CBCG, VectorAlgebra.scale(CGCD,-1) );//Bisector of the HB1-CB-HB2 angle
        double GPerpendicular[] = VectorAlgebra.cross( CBCG, CGCD );//Perpendicular to the CA-CB-CG plane
        double DBisector[] = VectorAlgebra.average( CGCD, VectorAlgebra.scale(CDN,-1) );//Bisector of the HB1-CB-HB2 angle
        double DPerpendicular[] = VectorAlgebra.cross( CGCD, CDN );//Perpendicular to the CA-CB-CG plane

        //We scale these vectors so they can be used as components of the C-H bond vectors
        BBisector = VectorAlgebra.scale( BBisector, HBDist * (double)Math.cos(HBAngle/2) / VectorAlgebra.norm(BBisector) );
        GBisector = VectorAlgebra.scale( GBisector, HGDist* (double)Math.cos(HGAngle/2) / VectorAlgebra.norm(GBisector) );
        DBisector = VectorAlgebra.scale( DBisector, HDDist* (double)Math.cos(HDAngle/2) / VectorAlgebra.norm(DBisector) );
        BPerpendicular = VectorAlgebra.scale( BPerpendicular, HBDist * (double)Math.sin(HBAngle/2) / VectorAlgebra.norm(BPerpendicular) );
        GPerpendicular = VectorAlgebra.scale( GPerpendicular, HGDist * (double)Math.sin(HGAngle/2) / VectorAlgebra.norm(GPerpendicular) );
        DPerpendicular = VectorAlgebra.scale( DPerpendicular, HDDist * (double)Math.sin(HDAngle/2) / VectorAlgebra.norm(DPerpendicular) );


        //Calculate the hydrogen coordinates
        //Using the AMBER hydrogen-naming system (HB1, HB2 etc. as in the AA templates)
        double HB1Coord[] = VectorAlgebra.add( CBCoord, VectorAlgebra.add(BBisector, BPerpendicular) );
        double HB2Coord[] = VectorAlgebra.add( CBCoord, VectorAlgebra.subtract(BBisector, BPerpendicular) );
        double HG1Coord[] = VectorAlgebra.add( CGCoord, VectorAlgebra.add(GBisector, GPerpendicular) );
        double HG2Coord[] = VectorAlgebra.add( CGCoord, VectorAlgebra.subtract(GBisector, GPerpendicular) );
        double HD1Coord[] = VectorAlgebra.add( CDCoord, VectorAlgebra.subtract(DBisector, DPerpendicular) );
        double HD2Coord[] = VectorAlgebra.add( CDCoord, VectorAlgebra.add(DBisector, DPerpendicular) );


        if ( res.getAtomIndexByName("HB1") != -1 ){//HB1-HB2 naming

            System.arraycopy( HB1Coord, 0, res.coords, 3*res.getAtomIndexByName("HB1"), 3 );
            System.arraycopy( HB2Coord, 0, res.coords, 3*res.getAtomIndexByName("HB2"), 3 );
            System.arraycopy( HG1Coord, 0, res.coords, 3*res.getAtomIndexByName("HG1"), 3 );
            System.arraycopy( HG2Coord, 0, res.coords, 3*res.getAtomIndexByName("HG2"), 3 );
            System.arraycopy( HD1Coord, 0, res.coords, 3*res.getAtomIndexByName("HD1"), 3 );
            System.arraycopy( HD2Coord, 0, res.coords, 3*res.getAtomIndexByName("HD2"), 3 );

        }
        else if ( res.getAtomIndexByName("HB2") != -1 ){//HB2-HB3 naming

            System.arraycopy( HB1Coord, 0, res.coords, 3*res.getAtomIndexByName("HB2"), 3 );
            System.arraycopy( HB2Coord, 0, res.coords, 3*res.getAtomIndexByName("HB3"), 3 );
            System.arraycopy( HG1Coord, 0, res.coords, 3*res.getAtomIndexByName("HG2"), 3 );
            System.arraycopy( HG2Coord, 0, res.coords, 3*res.getAtomIndexByName("HG3"), 3 );
            System.arraycopy( HD1Coord, 0, res.coords, 3*res.getAtomIndexByName("HD2"), 3 );
            System.arraycopy( HD2Coord, 0, res.coords, 3*res.getAtomIndexByName("HD3"), 3 );
        }
        else if ( res.getAtomIndexByName("2HB") != -1 ){//2HB-3HB naming

            System.arraycopy( HB1Coord, 0, res.coords, 3*res.getAtomIndexByName("2HB"), 3 );
            System.arraycopy( HB2Coord, 0, res.coords, 3*res.getAtomIndexByName("3HB"), 3 );
            System.arraycopy( HG1Coord, 0, res.coords, 3*res.getAtomIndexByName("2HG"), 3 );
            System.arraycopy( HG2Coord, 0, res.coords, 3*res.getAtomIndexByName("3HG"), 3 );
            System.arraycopy( HD1Coord, 0, res.coords, 3*res.getAtomIndexByName("2HD"), 3 );
            System.arraycopy( HD2Coord, 0, res.coords, 3*res.getAtomIndexByName("3HD"), 3 );
        }
        else{
            throw new Error("Can't parse atom names for " + res.fullName);
        }

    }
    
    
    
    
}
