/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author mhall44
 */
public class HardCodedResidueInfo {
    /* This version of OSPREY aims to make residue-type-dependent operations as generic as possible,
     * appropriate for D- and L-amino acids (natural or not) as well as non-amino acid residues
     * However, we sometimes require information that is not in the template files
     * The idea of this class is to get that stuff all in one place
     * to facilitate adding support for more residue types in the future
     * Also we may want to transfer information from here into some kind of template files in the future
     */
    
    
    public static String[] possibleBBAtoms = new String[] {
        "N", "H", "CA", "C", "O", "OXT", "H1", "H2", "H3"
    };
    //BB atoms, which should stay still in mutations and should be moved in perturbations.
    //We'll move HA with the sidechain, so it's not included here.  
    
    public static LinkedHashMap<String,String> three2one = null;

    static {
        initThree2One();
    }
    
    
    public static void initThree2One(){
            three2one = new LinkedHashMap<String,String>();
            three2one.put("ALA","A");
            three2one.put("CYS","C");
            three2one.put("ASP","D");
            three2one.put("GLU","E");
            three2one.put("PHE","F");
            three2one.put("GLY","G");
            three2one.put("HIS","H");
            three2one.put("ILE","I");
            three2one.put("LYS","K");
            three2one.put("LEU","L");
            three2one.put("MET","M");
            three2one.put("ASN","N");
            three2one.put("PRO","P");
            three2one.put("GLN","Q");
            three2one.put("ARG","R");
            three2one.put("SER","S");
            three2one.put("THR","T");
            three2one.put("VAL","V");
            three2one.put("TRP","W");
            three2one.put("TYR","Y");
    }

    public static String getOneLet(String aa3Name){
            String res = three2one.get(aa3Name);
            if (res == null)
                    res = "X";
            return res;
    }
        
        
        
    //Here's some stuff we need to mutate amino acid protein residues
    public static boolean canMutateTo(ResidueTemplate templ){
        //do we currently support mutations to the given amino-acid type?
        if(templ.templateRes.coords==null)
            return false;
        if(!hasAminoAcidBB(templ.templateRes))
            return false;
        
        return true;//can currently mutate to any amino acid (D or L, naturally occurring sidechain or not) 
        //whose sidechain attaches only to CA and for which we have template coords
    }
    
    
    public static ArrayList<String> listBBAtomsForMut(ResidueTemplate newTemplate, ResidueTemplate oldTemplate){
        //list the backbone atom names shared between the residue templates
        //(in the sense of atoms that don't move at all when we mutate--
        //basically backbone atoms that are present both before and after the mutation)
        
        if(!canMutateTo(newTemplate))
            throw new UnsupportedOperationException("ERROR: Can't currently mutate to "+newTemplate.name);
                        
        ArrayList<String> ans = new ArrayList<>();
        
        for(String atName : possibleBBAtoms){
            if(newTemplate.templateRes.getAtomIndexByName(atName) != -1){//atom name appears in first template
                if(oldTemplate.templateRes.getAtomIndexByName(atName) != -1){
                    ans.add(atName);
                }
            }
        }
        
        return ans;
    }
        
        
    public static int[][] findMutAlignmentAtoms(ResidueTemplate template1, ResidueTemplate template2){
        /*List indices of atoms in the two residue templates that can be aligned to perform a mutation
         * indices in first template go in ans[0], for second template go in ans[1]
         * 
         * Starting off with amino acids (L- or D- both supported, natural sidechains or not),
         * add more later (until then, trying to mutate a non-amino acid residue will cause an error though)
         * 
         * The alignment will use RigidBodyMotion.superimposingMotion, which prioritizes alignment of earlier atoms
         * so order matters
         */
        
        int mutAtoms1[] = getTemplateMutAtoms(template1);
        int mutAtoms2[] = getTemplateMutAtoms(template2);
        
        if(template1.name.equalsIgnoreCase("PRO") || template2.name.equalsIgnoreCase("PRO")){
            mutAtoms1 = getTemplateProMutAtoms(template1);
            mutAtoms2 = getTemplateProMutAtoms(template2);
        }
                
        for(int ma[] : new int[][] {mutAtoms1,mutAtoms2} ){
            for( int atNum : ma ){
                if(atNum==-1){//some atom(s) not found
                    throw new UnsupportedOperationException("ERROR: Mutation from "
                            +template1.name+" to "+template2.name+" not supported yet");
                }
            }
        }
        
        return new int[][] {mutAtoms1,mutAtoms2};
    }
            
    private static int[] getTemplateMutAtoms(ResidueTemplate template){
        //Get atoms to use for alignment
        int N = template.templateRes.getAtomIndexByName("N");
        int CA = template.templateRes.getAtomIndexByName("CA");
        int C = template.templateRes.getAtomIndexByName("C");
        
        return new int[] {CA,N,C};
    }
    
    //If mutating to or from PRO need to make sure CD is aligned to H (its non-Pro equivalent)
    //C and O will be OK...they are copied exactly
    private static int[] getTemplateProMutAtoms(ResidueTemplate template){
        //Get atoms to use for alignment
        int N = template.templateRes.getAtomIndexByName("N");
        int CA = template.templateRes.getAtomIndexByName("CA");
        
        if(template.name.equalsIgnoreCase("PRO")){//the PRO itself...use CD
            int CD = template.templateRes.getAtomIndexByName("CD");
            return new int[] {CA,N,CD};
        }
        else {//get H to align to CD
            int H = template.templateRes.getAtomIndexByName("H");
            return new int[] {CA,N,H};
        }
    }
    
    
    //The following functions are intended to mark inter-residue bonds between amino acids
    //We're assuming for now that non-amino acid residues aren't bonded to anything,
    //but could change that by changing just these functions
    //we might want a template file to tell this stuff eventually

    static double maxBondDist = 3;//liberal upper-bound on inter-residue bond length.  Given in angstroms

    
    //this version marks all the inter-residue bonds in a molecule
    public static void markInterResBonds(Molecule molec){
        
        
        //first, peptide bonds.  Must be between consecutive residues
        for(int resNum=0; resNum<molec.residues.size()-1; resNum++){
            Residue res1 = molec.residues.get(resNum);
            Residue res2 = molec.residues.get(resNum+1);
            tryToMakePeptideBond(res1,res2);
        }
        
        //there could also be disulfide bonds.  Look for any cysteins that are close enough
        for(Residue res1 : molec.residues){
            if(res1.fullName.startsWith("CYX")){//unprotonated cysteine...could be in a disulfide bond
                
                for(Residue res2 : molec.residues){
                    if(res2.indexInMolecule<res1.indexInMolecule){//don't double-count
                        if(res2.fullName.startsWith("CYX")){
                            tryToMakeDisulfideBond(res1,res2);
                        }
                    }
                }
            }
        }
        
        //now all inter-res bonds are assumed to be made
        for(Residue res : molec.residues)
            res.interResBondsMarked = true;
    }
    
    //this version just tries to make inter-res bonds involving the specified residue
    //it's intended for use right after res is mutated, to reconnect it to the rest of the molecule
    public static void reconnectInterResBonds(Residue res){
        
        //first try to make peptide bonds
        if(res.indexInMolecule>0){//res is not the first residue
            Residue prevRes = res.molec.residues.get(res.indexInMolecule-1);
            tryToMakePeptideBond(prevRes,res);
        }
        if(res.indexInMolecule<res.molec.residues.size()-1){//res is not the last residue
            Residue nextRes = res.molec.residues.get(res.indexInMolecule+1);
            tryToMakePeptideBond(res,nextRes);
        }
        
        
        //try to make any disulfide bonds, if res is an unprotonated cysteine
        if(res.fullName.startsWith("CYX")){
                
            for(Residue res2 : res.molec.residues){
                if(res2!=res){
                    if(res2.fullName.startsWith("CYX")){
                        tryToMakeDisulfideBond(res,res2);
                    }
                }
            }
        }
        
        //OK that should be all for now
        res.interResBondsMarked = true;
    }
    
    public static void tryToMakePeptideBond(Residue res1, Residue res2){
        //Given consecutive residues res1 and res2, make a peptide bond between them if appropriate
        if(res1.template==null || res2.template==null) {
            throw new RuntimeException("ERROR: Trying to peptide-bond residue without template");
        }
        if( (hasAminoAcidBB(res1) && hasAminoAcidBB(res2)) ){//can only peptide-bond amino acids
            
            int CIndex = res1.getAtomIndexByName("C");
            int NIndex = res2.getAtomIndexByName("N");
            
            //Get distance between these atoms
            double CNDist = VectorAlgebra.distance(res1.coords, CIndex, res2.coords, NIndex);
            if(CNDist<maxBondDist){
                Atom C = res1.atoms.get(CIndex);
                Atom N = res2.atoms.get(NIndex);
                C.addBond(N);
            }
        }
    }
    
    
    public static void tryToMakeDisulfideBond(Residue res1, Residue res2){
        //Given CYX residues res1 and res2, make a disulfide bond between them if appropriate
        if(res1.template==null || res2.template==null)
            throw new RuntimeException("ERROR: Trying to disulfide-bond residue without template");
        
        int SIndex1 = res1.getAtomIndexByName("SG");
        int SIndex2 = res2.getAtomIndexByName("SG");
        
        if(SIndex1==-1 || SIndex2==-1)
            throw new RuntimeException("ERROR: Trying to disulfide-bond residue without SG");

        //Get distance between these atoms
        double SSDist = VectorAlgebra.distance(res1.coords, SIndex1, res2.coords, SIndex2);
        if(SSDist<maxBondDist){
            Atom S1 = res1.atoms.get(SIndex1);
            Atom S2 = res2.atoms.get(SIndex2);
            S1.addBond(S2);
        }
    }
    
    
    public static boolean hasAminoAcidBB(Residue res){
        //does res have the three main amino-acid backbone atoms (N,CA, and C)?
        //This would allow it to peptide-bond
        if( res.getAtomIndexByName("N")>=0 && res.getAtomIndexByName("CA")>=0 
                && res.getAtomIndexByName("C")>=0){
            return true;//found them all
        }
        return false;//didn't find them all
    }
    
    
    public static String getTemplateName(Residue res){
        //some residues may have template names different than what's in the full name for the residue
        //we identify those here
        //we cover HIS and CYS; we assume templates follow the HID/HIE/HIP and CYS/CYX conventions as in
        //all_amino94.in
        
        String resName = res.fullName.substring(0,3);
        
        if( resName.equalsIgnoreCase("HIS") ){

            int HDCount = 0;//Count delta and epsilon hydrogens
            int HECount = 0;

            for( Atom at : res.atoms ){
                if( at.name.contains("HD") )
                    HDCount++;
                else if(at.name.contains("HE"))
                    HECount++;
            }
            
            //Determine protonation state, assign template name accordingly
            if( ( HDCount == 1 ) && ( HECount == 2 ) )
                return "HIE";
            else if( ( HDCount == 2 ) && ( HECount == 1 ) )
                return "HID";
            else if( ( HDCount == 2 ) && ( HECount == 2 ) )
                return "HIP";
            else{
                throw new RuntimeException("ERROR: Invalid protonation state for " + res.fullName );
            }
        }
        else if( resName.equalsIgnoreCase("CYS") ){
            
            if(res.getAtomIndexByName("HG")>=0)//there is a gamma hydrogen
                return "CYS";
            else
                return "CYX";
        }
        else//by default, the template name is just the first three letters of the full name
            return resName;
    }
    
    
    //Stuff about protonation-dependent rotamers for K* (Hie,Hid,etc.) could go in this class too...
    
}
