package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ResidueInteractions.Pair;

/** Residue interactions generator */
public class ResInterGen {
	
	private final SimpleConfSpace confSpace;
	private final ResidueInteractions inters;
	
	public static ResInterGen of(SimpleConfSpace confSpace) {
		return new ResInterGen(confSpace);
	}
	
	private ResInterGen(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
		this.inters = new ResidueInteractions();
	}
	
	public ResidueInteractions make() {
		return inters;
	}
	
	public String getResNum(int pos) {
		return confSpace.positions.get(pos).resNum;
	}
	
	public ResInterGen addIntra(int pos) {
		return addIntra(pos, Pair.IdentityWeight);
	}
	
	public ResInterGen addIntra(int pos, double weight) {	
		inters.addSingle(getResNum(pos), weight);
		return this;
	}
	
	public ResInterGen addIntras(RCTuple frag) {
		return addIntras(frag, Pair.IdentityWeight);
	}
	
	public ResInterGen addIntras(RCTuple frag, double weight) {
		for (int i=0; i<frag.size(); i++) {
			int pos = frag.pos.get(i);
			inters.addSingle(getResNum(pos), weight);
		}
		return this;
	}
	
	public ResInterGen addInter(int pos1, int pos2) {
		return addInter(pos1, pos2, Pair.IdentityWeight);
	}
	
	public ResInterGen addInter(int pos1, int pos2, double weight) {
		inters.addPair(getResNum(pos1), getResNum(pos2), weight);
		return this;
	}
	
	public ResInterGen addInters(RCTuple frag) {
		return addInters(frag, Pair.IdentityWeight);
	}
	
	public ResInterGen addInters(RCTuple frag, double weight) {
		for (int i=0; i<frag.size(); i++) {
			String resNum1 = getResNum(frag.pos.get(i));
			for (int j=0; j<i; j++) {
				inters.addPair(resNum1, getResNum(frag.pos.get(j)), weight);
			}
		}
		return this;
	}
	
	public ResInterGen addShell(int pos) {
		return addShell(pos, Pair.IdentityWeight);
	}
	
	public ResInterGen addShell(int pos, double weight) {
		String resNum = getResNum(pos);
		for (String resNumShell : confSpace.shellResNumbers) {
			inters.addPair(resNum, resNumShell, weight);
		}
		return this;
	}
	
	public ResInterGen addShell(RCTuple frag) {
		return addShell(frag, Pair.IdentityWeight);
	}
	
	public ResInterGen addShell(RCTuple frag, double weight) {
		for (int i=0; i<frag.size(); i++) {
			String resNum = getResNum(frag.pos.get(i));
			for (String resNumShell : confSpace.shellResNumbers) {
				inters.addPair(resNum, resNumShell, weight);
			}
		}
		return this;
	}
	
	public ResInterGen add(double offset) {
		inters.addOffset(offset);
		return this;
	}
}
