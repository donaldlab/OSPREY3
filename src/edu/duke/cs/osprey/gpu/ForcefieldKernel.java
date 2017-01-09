package edu.duke.cs.osprey.gpu;

import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;

public interface ForcefieldKernel {
	
	BigForcefieldEnergy getForcefield();
	BigForcefieldEnergy.Subset getSubset();
	boolean setSubset(BigForcefieldEnergy.Subset subset);
	void uploadCoordsAsync();
	void runAsync();
	void waitForGpu();
	double downloadEnergySync();
	void cleanup();
}
