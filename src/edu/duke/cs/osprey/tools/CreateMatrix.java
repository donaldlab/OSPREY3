package edu.duke.cs.osprey.tools;

/**
 * Creates matrices needed for different purposes in OSPREY
 * @author pablo
 *
 */
public class CreateMatrix {
	// Creates a 4d rotamer matrix of the form [resI][rotR][resJ][rotS] = initializationValue
	public static double [][][][] create4DRotMatrix(int numRes, int rotsPerRes[], double initializeValue){
		double [][][][] rot4DMat = new double[numRes][][][];
		for (int resI = 0; resI < numRes; resI++){
			rot4DMat[resI] = new double[rotsPerRes[resI]][][];
			for (int rotIR = 0; rotIR < rotsPerRes[resI]; rotIR++){
				rot4DMat[resI][rotIR] = new double[numRes][]; 
				for(int resJ = 0; resJ < numRes; resJ++){
					rot4DMat[resI][rotIR][resJ] = new double [rotsPerRes[resJ]];
					for(int rotJS = 0; rotJS < rotsPerRes[resJ]; rotJS++){
						rot4DMat[resI][rotIR][resJ][rotJS] = initializeValue;
					}					
				}
			}
		}
		return rot4DMat;
	}
	// Creates a 2d rotamer matrix of the form [resI][rotR] = initializationValue
	public static double [][] create2DRotMatrix(int numRes, int rotsPerRes[], double initializeValue){
		double [][] rot2DMat = new double[numRes][];
		for (int resI = 0; resI < numRes; resI++){
			rot2DMat[resI] = new double[rotsPerRes[resI]];
			for (int rotIR = 0; rotIR < rotsPerRes[resI]; rotIR++){
				rot2DMat[resI][rotIR] = initializeValue; 
				
			}
		}
		return rot2DMat;
	}
	// Creates a 3d rotamer matrix of the form [senderRes][recvRes][recvRot]
	public static double [][][] create3DMsgMat(int numRes, int rotsPerRes[], double initializeValue){
		double [][][] msg3DMat = new double[numRes][][];
		for (int sendRes = 0; sendRes < numRes; sendRes ++){
			msg3DMat[sendRes] = new double[numRes][];
			for (int recvRes = 0; recvRes < numRes; recvRes++){
				msg3DMat[sendRes][recvRes] = new double[rotsPerRes[recvRes]];
				for (int recvRot = 0; recvRot < rotsPerRes[recvRes]; recvRot++){
					msg3DMat[sendRes][recvRes][recvRot] = initializeValue;
				}
			}
		}
		return msg3DMat;
	}

}