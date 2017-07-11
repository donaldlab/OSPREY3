package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;

/**
 * Maintains the design positions and residue conformations for a list of strands.
 * 
 * Also creates molecules in conformations on-demand.
 */
public class SimpleConfSpace implements Serializable {
	
	public static class Builder {
		
		private SimpleConfSpace confSpace;
		
		/**
		 * A residue is included in the steric shell if any of its atoms lies within
		 * shellDist angstroms of any atom in any flexible residue.
		 */
		private double shellDist = Double.POSITIVE_INFINITY;
		
		public Builder() {
			confSpace = new SimpleConfSpace();
		}
		
		public Builder addStrand(Strand strand, StrandFlex ... flexType) {
			addStrand(strand, Arrays.asList(flexType));
			return this;
		}
		
		public Builder addStrand(Strand strand, List<StrandFlex> flexTypes) {
			confSpace.addStrand(strand, flexTypes);
			return this;
		}
		
		public Builder addStrands(Strand ... strands) {
			for (Strand strand : strands) {
				addStrand(strand);
			}
			return this;
		}
		
		public Builder setShellDistance(double val) {
			shellDist = val;
			return this;
		}
		
		public SimpleConfSpace build() {
			
			// make sure we have some design positions
			if (confSpace.positions.isEmpty()) {
				throw new IllegalStateException("ConfSpace has no design positions, try adding some strand flexibility");
			}
			
			confSpace.makeShell(shellDist);
			confSpace.countResConfs();
			return confSpace;
		}
	}
	
	public static class Position implements Serializable {
		
		public final int index;
		public final Strand strand;
		public final String resNum;
		public final Strand.ResidueFlex resFlex;
		public final List<ResidueConf> resConfs;
		
		public Position(int index, Strand strand, Residue res) {
			this.index = index;
			this.strand = strand;
			this.resNum = res.getPDBResNumber();
			this.resFlex = strand.flexibility.get(resNum);
			this.resConfs = new ArrayList<>();
		}
	}
	
	public static class ResidueConf implements Serializable {
		
		public static enum Type {
			
			Library('L'),
			WildType('W');
			
			public final char letter;
			
			private Type(char letter) {
				this.letter = letter;
			}
			
			public static Type getFromLetter(char letter) {
				for (Type type : values()) {
					if (type.letter == letter) {
						return type;
					}
				}
				return null;
			}
		}
		
		public static interface PostTemplateModifier extends Serializable {
			void modify(Residue res);
		}
		
		public final int index;
		public final ResidueTemplate template;
		public final Type type;
		public final Integer rotamerIndex;
		
		public PostTemplateModifier postTemplateModifier;
                
                public HashMap<String,double[]> dofBounds;//for DEEPer need this, else it's just a rotamer
                //for EPIC it is necessary to have some kind of persistent representation of the DOFs.  Using names here
                //so can dynamically create dihedral DOFs, etc. without getting confused
		                
		public ResidueConf(Position pos, int index, ResidueTemplate template, Type type, HashMap<String,double[]> bbVoxel) {
			this(pos, index, template, type, null, bbVoxel);
		}
		
		public ResidueConf(Position pos, int index, ResidueTemplate template, Type type, Integer rotamerIndex, HashMap<String,double[]> bbVoxel) {
			
			if (template.getNumRotamers() > 0 && rotamerIndex == null) {
				throw new IllegalArgumentException("template " + template.name + " has rotamers, so rotamer index is required");
			}
			
			this.index = index;
			this.template = template;
			this.type = type;
			this.rotamerIndex = rotamerIndex;
			
			this.postTemplateModifier = null;
                        
                        dofBounds = sidechainDOFBounds(pos,template,rotamerIndex);
                        dofBounds.putAll(bbVoxel);
		}

		public void updateResidue(GenericResidueTemplateLibrary templateLib, Residue res) {
			
			// HACKHACK: make sure prolines have puckers
			if (res.template.name.equalsIgnoreCase("PRO") || template.name.equalsIgnoreCase("PRO")) {
				res.pucker = new ProlinePucker(templateLib, res);
			}
			
			if (postTemplateModifier != null) {
				postTemplateModifier.modify(res);
			}
			
			ResidueTypeDOF.switchToTemplate(templateLib, res, template, false);
		}
		
		public String getRotamerCode() {
			StringBuilder buf = new StringBuilder();
			buf.append(type.letter);
			if (rotamerIndex != null) {
				buf.append(rotamerIndex);
			}
			return buf.toString();
		}
		
		@Override
		public String toString() {
			return template.name + " " + getRotamerCode();
		}
                
                public boolean isParametricallyIncompatibleWith(ResidueConf rc2){
                    //Two RCs are parametrically incompatible if and only if there is a DOF that they share
                    //for which they have different intervals
                    //Thus, a conformation has a well-defined voxel if and only if it contains no parametrically
                    //incompatible pairs of RCs

                    final double tol = 1e-8;
                    for(String dofName1 : dofBounds.keySet()){
                        for(String dofName2 : rc2.dofBounds.keySet()){
                            if(dofName1.equalsIgnoreCase(dofName2)){
                                double bounds1[] = dofBounds.get(dofName1);
                                double bounds2[] = dofBounds.get(dofName2);
                                for(int a=0; a<2; a++){//look for min or max not matching
                                    if( Math.abs( bounds1[a]-bounds2[a] ) > tol )
                                        return true;
                                }
                            }
                        }
                    }

                    //no incompatibility found!
                    return false;
                }
	}
	
	/**
	 * What kind of degree of freedom?
	 * 
	 * Used mainly to decide if we can use GPU CCD implementations or not,
	 * which only support dihedral angles for now.
	 */
	public static enum DofTypes {
		
		// NOTE: order is important for combine()
		None,
		OnlyDihedrals,
		Any;
		
		public static DofTypes combine(DofTypes a, DofTypes b) {
			if (a == b) {
				return a;
			} else if (a == null && b != null) {
				return b;
			} else if (a != null && b == null) {
				return a;
			} else {
				if (a.ordinal() > b.ordinal()) {
					return a;
				} else {
					return b;
				}
			}
		}
	}
	
	
	/** The strands */
	public final List<Strand> strands;
	
	/** The design positions */
	public final List<Position> positions;
	
	/** The residue numbers of the shell residues */
	public final Set<String> shellResNumbers;
	
	/** The flexibility of each strand */
	public final Map<Strand,List<StrandFlex>> strandFlex; // yeah, map on instance, not identity
	
	private Map<Strand,Set<String>> shellResNumbersByStrand;
	private int[] numResConfsByPos;
        
	private SimpleConfSpace() {
		strands = new ArrayList<>();
		positions = new ArrayList<>();
		shellResNumbers = new HashSet<>();
		strandFlex = new HashMap<>();
		
		shellResNumbersByStrand = new HashMap<>();
		numResConfsByPos = null;
	}
	
	private void addStrand(Strand strand, List<StrandFlex> flexTypes) {
		
		// add the strand
		strands.add(strand);
		strandFlex.put(strand, flexTypes);
		
		// make the positions
		for (String resNum : strand.flexibility.getFlexibleResidueNumbers()) {
			Residue res = strand.mol.getResByPDBResNumber(resNum);
			Strand.ResidueFlex resFlex = strand.flexibility.get(resNum);
			
			// make the pos
			Position pos = new Position(positions.size(), strand, res);
			positions.add(pos);
			
			// make residue confs from library rotamers
			for (String resType : resFlex.resTypes) {
				makeResidueConfsFromTemplate(pos, strand.templateLib.getTemplate(resType), ResidueConf.Type.Library);
			}
			
			// make residue confs from wild type rotamers
			if (resFlex.addWildTypeRotamers) {
				
				List<Residue> residues = new ArrayList<>();
				residues.add(res);
				residues.addAll(strand.mol.getAlternates(res.indexInMolecule));
				
				makeResidueConfsFromTemplate(pos, ResidueTemplate.makeFromResidueConfs(residues), ResidueConf.Type.WildType);
			}
		}
	}
	
	private void makeResidueConfsFromTemplate(Position pos, ResidueTemplate template, ResidueConf.Type type) {
		
            //make one RC for each (backbone voxel, rotamer) pair
            for(HashMap<String,double[]> bbState : listBackboneVoxels(pos)){ 
            
		if (template.name.equalsIgnoreCase("PRO")) {
			
			// HACKHACK: add one conf for each proline pucker
			for (ProlinePucker.Direction dir : ProlinePucker.Direction.values()) {
				ResidueConf resConf = new ResidueConf(
                                        pos,
					pos.resConfs.size(),
					template,
					type,
					dir.ordinal(),
                                        bbState
				);
				resConf.postTemplateModifier = (res) -> res.pucker.apply(dir);
				pos.resConfs.add(resConf);
			}
			
		} else if (template.getNumRotamers() <= 0) {
			
			// make one template for the library template
			pos.resConfs.add(new ResidueConf(
                                pos,
				pos.resConfs.size(),
				template,
				type,
                                bbState
			));
			
		} else {
			
			// make a template for each rotamer
			for (int rotamerIndex=0; rotamerIndex<template.getNumRotamers(); rotamerIndex++) {
				pos.resConfs.add(new ResidueConf(
                                        pos,
					pos.resConfs.size(),
					template,
					type,
					rotamerIndex,
                                        bbState
				));
			}
		}
            }
	}
	
	private void makeShell(double shellDist) {
		
		shellResNumbers.clear();
		shellResNumbersByStrand.clear();
	
		for (Strand strand : strands) {
			
			Set<String> resNumbers = new HashSet<>();
			
			for (String staticResNum : strand.flexibility.getStaticResidueNumbers()) {
				Residue staticRes = strand.mol.getResByPDBResNumber(staticResNum);
				
				for (String flexResNum : strand.flexibility.getFlexibleResidueNumbers()) {
					Residue flexRes = strand.mol.getResByPDBResNumber(flexResNum);
					
					if (staticRes.distanceTo(flexRes) <= shellDist) {
						resNumbers.add(staticResNum);
						break;
					}
				}
			}
			
			shellResNumbers.addAll(resNumbers);
			shellResNumbersByStrand.put(strand, resNumbers);
		}
	}
	
	private void countResConfs() {
		numResConfsByPos = new int[positions.size()];
		for (int i=0; i<positions.size(); i++) {
			numResConfsByPos[i] = positions.get(i).resConfs.size();
		}
	}
	
	/** Gets the number of residue confs per position */
	public int[] getNumResConfsByPos() {
		return numResConfsByPos;
	}
	
	/** Gets the number of residue confs for a position */
	public int getNumResConfs(int pos) {
		return positions.get(pos).resConfs.size();
	}
	
	/** Gets the total number of residue confs for all positions */
	public int getNumResConfs() {
		int count = 0;
		for (int pos=0; pos<positions.size(); pos++) {
			count += numResConfsByPos[pos];
		}
		return count;
	}
	
	/** Gets the total number of residue conf pairs for all positions */
	public int getNumResConfPairs() {
		int count = 0;
		for (int pos1=0; pos1<positions.size(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				for (int rc1=0; rc1<numResConfsByPos[pos2]; rc1++) {
					count += numResConfsByPos[pos1];
				}
			}
		}
		return count;
	}
	
	/** @see #makeMolecule(RCTuple) */
	/*public ParametricMolecule makeMolecule(ScoredConf conf) {
		return makeMolecule(conf.getAssignments());
	}*/
        public BoundedParametricMolecule makeBoundedParametricMolecule(ScoredConf conf) {
		return makeBoundedParametricMolecule(conf.getAssignments());
	}
        //With DEEPer and CATS it is much less messy if we can generate the ParametricMolecule and its bounds together
	
	/** @see #makeMolecule(RCTuple) */
	/*public ParametricMolecule makeMolecule(int[] conf) {
		return makeMolecule(new RCTuple(conf));
	}*/
        public BoundedParametricMolecule makeBoundedParametricMolecule(int[] conf) {
		return makeBoundedParametricMolecule(new RCTuple(conf));
	}
	
	/**
	 * create a new {@link ParametricMolecule} in the specified conformation
	 * for analysis (e.g., energy calculation, minimization)
	 * 
	 * To increase stability of analysis, each analysis should be conducted
	 * with a new molecule instance. this completely prevents roundoff error
	 * from accumulating across separate analyses. 
	 */
	public BoundedParametricMolecule makeBoundedParametricMolecule(RCTuple conf) {
		
		// make the molecule from the strands (ignore alternates)
		Molecule mol = new Molecule();
		for (Strand strand : strands) {
			for (Residue res : strand.mol.residues) {
				res = new Residue(res);
				res.molec = mol;
				res.indexInMolecule = mol.residues.size();
				mol.residues.add(res);
			}
		}
		HardCodedResidueInfo.markInterResBonds(mol);
		
		// mutate to the conf templates, and figure out what conformational DOFs are specified by the conf
                HashSet<String> confDOFNames = new HashSet<>();//names of DOFs specified by the conf
		for (int i=0; i<conf.size(); i++) {
			
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			Residue res = mol.getResByPDBResNumber(pos.resNum);
			
			resConf.updateResidue(pos.strand.templateLib, res);
                        //since we always switch to a new template before starting each minimization,
                        //no need to standardize mutatable res at the beginning of the design
                        
                        for(String dofName : resConf.dofBounds.keySet())
                            confDOFNames.add(dofName);
		}
                
                
                //OK now apply all DOF vals including puckers
		
		// make all the DOFs
		List<DegreeOfFreedom> dofs = new ArrayList<>();
		
		// first, backbone flexibility: strand DOFs and DEEPer/CATS DOFs
		for (Strand strand : getConfStrands(conf)) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				//dofs.addAll(flex.makeDofs(strand));
                                for(DegreeOfFreedom dof : flex.makeDofs(strand, mol)){
                                    if(confDOFNames.contains(dof.getName())){
                                        //DEEPer and CATS DOFS may not involve all positions in a strand
                                        dofs.add(dof);
                                    }
                                }
			}
		}
		
		// then, residue conf DOFs
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			Residue res = mol.getResByPDBResNumber(pos.resNum);
			
			// make the residue DOFs
			Strand.ResidueFlex resFlex = pos.strand.flexibility.get(pos.resNum);
                        List<DegreeOfFreedom> contDihedralDOFs = resFlex.voxelShape.makeDihedralDOFs(res);
			dofs.addAll(contDihedralDOFs);
                        
                        //makeDihedralDOFs will not generate DOFs for a rigid rotamer
                        //so we apply them here instead
                        if(contDihedralDOFs.isEmpty()){
                            // pose the residue to match the rotamer
                            List<DegreeOfFreedom> dihedralDofs = new VoxelShape.Rect().makeDihedralDOFs(res);
                            for (int d=0; d<resConf.template.numDihedrals; d++)
                                    dihedralDofs.get(d).apply(resConf.template.getRotamericDihedrals(resConf.rotamerIndex, d));
                        }
                }
                
                //Figure out the bounds and apply the DOF values in the middle of the voxel
                DofBounds dofBounds = new DofBounds(dofs.size());//the bounds to use
                HashMap<String,Integer> name2Index = DegreeOfFreedom.nameToIndexMap(dofs);//map DOF names to their indices in DOFs
                HashSet<String> dofsAdded = new HashSet<>();
                
                for (int i=0; i<conf.size(); i++) {
                        Position pos = positions.get(conf.pos.get(i));
                        ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
                        HashMap<String,double[]> rcDOFBounds = resConf.dofBounds;
                        for(String DOFName : rcDOFBounds.keySet()){
                            int dofIndex = name2Index.get(DOFName);
                            double curDOFBounds[] = rcDOFBounds.get(DOFName);
                            if(dofsAdded.contains(DOFName)){//we have bounds on this DOF from another position...check consistency
                                if( Math.abs(dofBounds.getMax(dofIndex)-curDOFBounds[1])>1e-10 
                                        || Math.abs(dofBounds.getMin(dofIndex)-curDOFBounds[0])>1e-10 ){
                                    throw new RuntimeException("ERROR: Conformation has inconsistent DOF bounds"
                                            + " between different positions' RCs");
                                }
                            }
                            else {//record the bounds and apply the middle value
                                dofBounds.set(dofIndex, curDOFBounds[0], curDOFBounds[1]);
                                dofs.get(dofIndex).apply(0.5*(curDOFBounds[0]+curDOFBounds[1]));
                                dofsAdded.add(DOFName);
                            }
                        }
                }
		
		return new BoundedParametricMolecule( new ParametricMolecule(mol, dofs), dofBounds );
	}
	
	/** @see #isContinuouslyFlexible(RCTuple) */
	public boolean isContinuouslyFlexible(int[] conf) {
		return isContinuouslyFlexible(new RCTuple(conf));
	}
	
	/** Return true if the conformation has continuous degrees of freedom, false otherwise. */
	public boolean isContinuouslyFlexible(RCTuple conf) {
		for (int i=0; i<conf.size(); i++) {
                    Position pos = positions.get(conf.pos.get(i));
                    ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
                    for(double[] bounds : resConf.dofBounds.values()){
                        if(bounds[1]>bounds[0])
                            return true;
                    }
                }
                return false;
	}
	
	/** @see #makeBounds(RCTuple) */
	/*public DofBounds makeBounds(int[] conf, List<DegreeOfFreedom> dofs) {
		return makeBounds(new RCTuple(conf), dofs);
	}*/
	
	/** Calculate the bounds on all the degrees of freedom. */
        //MH: Since different tuples can have different CATS and DEEPer DOFs
        //we need to specify what DOFs we want to bound
	/*public DofBounds makeBounds(RCTuple conf, List<DegreeOfFreedom> dofs) {
            DofBounds ans;
            HashMap<String,Integer> name2Index = DegreeOfFreedom.nameToIndexMap(dofs);
            for (int i=0; i<conf.size(); i++) {
                    Position pos = positions.get(conf.pos.get(i));
                    ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
                    HashMap<String,double[]> rcDOFBounds = resConf.DOFBoundsMap;
                    for(String DOFName : rcDOFBounds){
                        ans[name2Index(DOFName)] = rcDOFBounds[DOFName];
                        //also check if there is already something in there by that name
                    }
            }
            return ans;
        }*/        
        //DEBUG!!!
        /*public DofBounds makeBounds(RCTuple conf) {
		
		// gather all the DOF bounds
		List<DofBounds> bounds = new ArrayList<>();
		
		// first, strand bounds
		for (Strand strand : getConfStrands(conf)) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				bounds.add(flex.makeBounds(strand));
			}
		}
		
		// then, residue conf bounds
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			if (resConf.rotamerIndex != null) {
				bounds.add(pos.resFlex.voxelShape.makeDihedralBounds(resConf.template, resConf.rotamerIndex));
			}
		}
		
		// TODO: make DOF bounds for non-residue continuous motions (eg, blocks, DEEPer)
		
		return DofBounds.concatenate(bounds);
	}*/
	
	private Set<Strand> getConfStrands(RCTuple conf) {
		Set<Strand> confStrands = new HashSet<>();
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			confStrands.add(pos.strand);
		}
		return confStrands;
	}

	public DofTypes getDofTypes() {
		
		DofTypes dofTypes = null;
		
		// check strand flex
		for (Strand strand : strands) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				dofTypes = DofTypes.combine(dofTypes, flex.getDofTypes());
			}
		}
		
		// check residue flex
		for (Position pos : positions) {
			dofTypes = DofTypes.combine(dofTypes, pos.resFlex.voxelShape.getDofTypes());
		}
		
		return dofTypes;
	}
        
        public int getNumPos(){
            return positions.size();
        }
        
        private static HashMap<String,double[]> sidechainDOFBounds(Position pos, ResidueTemplate template, Integer rotamerIndex){
            //get bounds on the sidechain dihedrals associated with a particular rotamer
            if(rotamerIndex==null)//no dihedrals
                return new HashMap<>();
            
            DofBounds bounds = pos.resFlex.voxelShape.makeDihedralBounds(template, rotamerIndex);
            Residue res = pos.strand.mol.getResByPDBResNumber(pos.resNum);

            HashMap<String,double[]> ans = new HashMap<>();
            for(int d=0; d<bounds.size(); d++){//in rigid case will not make bounds
                FreeDihedral dih = new FreeDihedral(res,d);
                ans.put(dih.getName(), new double[]{bounds.getMin(d),bounds.getMax(d)});
            }
            
            return ans;
        }
        
        
        private ArrayList<HashMap<String,double[]>> listBackboneVoxels(Position pos){
            ArrayList<HashMap<String,double[]> > ans = new ArrayList<>();
            
            for(StrandFlex flexType : strandFlex.get(pos.strand)){
                ArrayList<HashMap<String,double[]> > flexTypeVox = flexType.listBackboneVoxels(pos);
                if( (!flexTypeVox.isEmpty()) && (!ans.isEmpty()) ){
                    throw new RuntimeException("ERROR: Can't have multiple types of backbone flexibility for the same residue");
                    //not supported, because current backbone DOF implementations depend on the DOF block
                    //keeping track of the unperturbed backbone conformation, so if another type of motion
                    //changes that unperturbed bb conf, then there will be errors
                    //Mutations and sidechain dihedrals don't move the backbone so no issue there
                }
                else
                    ans = flexTypeVox;
            }
            
            if(ans.isEmpty())//no backbone flexibility
                return new ArrayList(Arrays.asList(new HashMap<>()));
            else
                return ans;
        }
}
