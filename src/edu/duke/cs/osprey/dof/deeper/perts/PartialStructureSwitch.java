/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper.perts;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.deeper.ResBBState;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;

import java.util.ArrayList;

/**
 *
 * This is a perturbation that substitutes in the backbone conformation from alternate structures
 * 
 * @author mhall44
 */
public class PartialStructureSwitch extends Perturbation {

    ArrayList<ArrayList<ResBBState>> altResConfs = new ArrayList<>();
    //residues (read from alternate PDB files) giving alternate backbone
    //conformations for each of the resDirectlyAffected
    //Indexed by structure (i.e., param val), then residue

    
    public PartialStructureSwitch(ArrayList<Residue> resDirectlyAffected, ArrayList<String> altConfPDBFiles) {
        super(resDirectlyAffected);
        
        altResConfs.add(null);//skip unperturbed structure
        
        //The alternate conformations are given by the residues with
        //the same numbers as resDirectlyAffected
        //in the altConfPDBFiles
        for(String altPDB : altConfPDBFiles){
            Molecule altMolec = new Strand.Builder(PDBIO.readFile(altPDB)).build().mol;
            ArrayList<ResBBState> altConf = new ArrayList<>();
            
            for(Residue origRes : resDirectlyAffected){
                Residue altRes = altMolec.getResByPDBResNumber(origRes.getPDBResNumber());
                altConf.add( new ResBBState(altRes) );
            }
            altResConfs.add(altConf);
        }
    }
    
    
    //for copying
	public PartialStructureSwitch(ArrayList<Residue> resDirectlyAffected, double curParamVal,
            int indexInBlock, PerturbationBlock block, ArrayList<ArrayList<ResBBState>> oldAltResConfs){
        super(resDirectlyAffected);
        this.curParamVal = curParamVal;
        this.indexInBlock = indexInBlock;
        this.block = block;
        
        for(ArrayList<ResBBState> arcOld : oldAltResConfs){
            ArrayList<ResBBState> arc = new ArrayList<>();
            if(arcOld == null) arc = null;
            else {
            	for(ResBBState rbs : arcOld)
            		arc.add(new ResBBState(rbs));
            }
            altResConfs.add(arc);
        }
    }
    
    
    @Override
    public boolean doPerturbationMotion(double paramVal) {
        
        int structNum = (int)paramVal;
        
        if(structNum>0){//0 indicates a state besides the original one
            for(int resNum=0; resNum<resDirectlyAffected.size(); resNum++){
                //just copy over backbone coordinates from other residue
                Residue res = resDirectlyAffected.get(resNum);
                
                ResBBState altConf = altResConfs.get(structNum).get(resNum);
                altConf.putInState(res);
                
                /*Residue altConf = altResConfs.get(structNum).get(resNum);
                
                int numAtoms = res.atoms.size();
                
                for(int atomNum=0; atomNum<numAtoms; atomNum++){
                    String atomName = res.atoms.get(atomNum).name;
                    if(isBBAtom(res,atomName)){//BB atom: copy coordinates
                        int atomNumAlt = altConf.getAtomIndexByName(atomName);
                        System.arraycopy(altConf.coords, 3*atomNumAlt, res.coords, 3*atomNum, 3);
                    }
                }    */            
            }
        }
        
        return true;//A partial structure switch will always put in structure as given.  
    }
    
    
    private static boolean isBBAtom(Residue res, String atomName){
        //Check if the atom should be moved as part of the backbone
        if(res.template.name.equalsIgnoreCase("PRO")){//CD moves with backbone
            if(atomName.equalsIgnoreCase("CD"))
                return true;
        }
        
        for(String BBName : HardCodedResidueInfo.possibleBBAtoms){
            if(atomName.equalsIgnoreCase(BBName))
                return true;
        }
        
        
        //doesn't match a known backbone atom name
        return false;
    }
    
    
    @Override
    public Perturbation copyForNewMolecule(Molecule mol, PerturbationBlock block){
        return new PartialStructureSwitch(
                Residue.equivalentInMolec(resDirectlyAffected, mol),
                curParamVal,
                indexInBlock,
                block,
                altResConfs
        );
    }
    
    
}
