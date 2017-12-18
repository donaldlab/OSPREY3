/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Set;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;

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
    
    public static Set<String> possibleBBAtomsLookup;
    
    //BB atoms, which should stay still in mutations and should be moved in perturbations.
    //We'll move HA with the sidechain, so it's not included here.  
    
    public static LinkedHashMap<String,String> three2one = null;
    public static LinkedHashMap<String,String> one2three = null;//reverse lookup

    static {
        initThree2One();
        
        one2three = new LinkedHashMap<String,String>();
        for(String threeLet : three2one.keySet())
            one2three.put(three2one.get(threeLet), threeLet);
        
        possibleBBAtomsLookup = new HashSet<>();
        for (String name : possibleBBAtoms) {
            possibleBBAtomsLookup.add(name);
        }
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




    public static double bondLengthUpperBound(int elementNum1, int elementNum2){
        //These are used to infer inter-res bonds, which cannot be specified directly by the template
        //because we don't know what other residues a given residue type will be bonded to
        //OSPREY 2 did this for all bonds; here we only use it for a few atoms
        //OSPREY 2 put in bonds if distance < 1.5 for S-H, 1.2 for other -H,
        //2.4 for S-S and any -Br, 2 for other -S, -Cl, -P,
        //otherwise 1.7

        //Moving 2 up to 2.2 in response to 2.01 A P-S bond in 1AK0 res 295

        if(Math.min(elementNum1, elementNum2)==1){
            //bond involving hydrogen.
            return Math.max(elementNum1,elementNum2)==16 ? 1.5 : 1.2;
        }

        switch(elementNum1){
            case 35:
                return 2.4;
            case 16:
                return (elementNum2==16 || elementNum2==35) ? 2.4 : 2.2;
            case 15: case 17:
                return (elementNum2==35) ? 2.4 : 2.2;
            default:
                switch(elementNum2){
                    case 35:
                        return 2.4;
                    case 15: case 16: case 17:
                        return 2.2;
                    default:
                        return 1.7;
                }
        }
    }
    

    public static InterResBondingTemplate inferInterResBonding(Residue res){
        //DEBUG!!  If no amino-acid backbone assuming no inter-res bonds
        if(HardCodedResidueInfo.hasAminoAcidBB(res)){
            if(res.fullName.startsWith("CYX"))
                return new InterResBondingTemplate.CysteineBondingTemplate();
            else
                return new InterResBondingTemplate.PeptideBondingTemplate();
        }
        else
            return new InterResBondingTemplate.NoBondingTemplate();
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
