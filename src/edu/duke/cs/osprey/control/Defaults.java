package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.FORCEFIELD;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.FileTools;

public class Defaults {
	
	/* IMPORTANT: rules for how to handle configuration without creating unmaintainable code:
	 * 
	 * For configuration to be effective and flexible, classes needing configuration
	 * should always allow the config to be passed in as arguments, set as class members, etc.
	 * Config for each object should only query these default values as a fallback
	 * when config is not provided at object creation time (via eg constructors, member setters).
	 * 
	 * It is crucially important that these default values remain perfectly constant
	 * over the course of each Osprey run.
	 * 
	 * Changing these default values during an Osprey run is strictly prohibited.
	 * 
	 * If a class instance requires configuration other than the default values,
	 * that class should support changing its internal configuration via constructors, member setters, etc.
	 * This is usually accomplished by allocating separate instances of the configuration objects
	 * (eg ForcefieldParams, ResidueTemplateLibrary) and assigning them to the class instance that needs config.
	 * 
	 * Under no circumstances should a class ONLY query configuration from the defaults.
	 * Such classes are not sufficiently flexible for use under different contexts.
	 * Classes should always allow configuration to be set by the caller of that class.
	 * 
	 * Once a class has been configured by a caller, that class should not think it "owns" that configuration.
	 * The caller "owns" that configuration. Classes should never modify their configuration objects.
	 * Callers could re-use configuration between classes, so classes must not break any assumptions
	 * callers may have made about the configuration objects.
	 *  
	 * If a class needs a configuration that is different from the one set by its caller,
	 * that class should allocate new configuration objects and configure them for internal use by the class.
	 * 
	 * GOOD EXAMPLE:
	 * 
	 * class GoodClass {
	 * 
	 *     public final Config config;
	 *     
	 *     public GoodClass() {
	 *         this(Defaults.config);
	 *     }
	 *     
	 *     public GoodClass(Config config) {
	 *         this.config = config;
	 *     }
	 *     
	 *     public void doTheThing() {
	 *         System.out.println(config.value);
	 *     }
	 * }
	 * 
	 * BAD EXAMPLE:
	 * 
	 * class BadClass {
	 *     
	 *     public void doTheThing() {
	 *         System.out.println(Defaults.config.value);
	 *     }
	 * }
	 */
	
	public static final FORCEFIELD forcefield;
	public static final ForcefieldParams forcefieldParams;
	public static final GenericResidueTemplateLibrary genericTemplateLibrary;
	
	static {
		
		forcefield = FORCEFIELD.AMBER;
		forcefieldParams = new ForcefieldParams();
		
		genericTemplateLibrary = new GenericResidueTemplateLibrary();
		genericTemplateLibrary.loadTemplateCoords(FileTools.readResource("/config/all_amino_coords.in"));
		genericTemplateLibrary.loadRotamerLibrary(FileTools.readResource("/config/LovellRotamer.dat"));
		genericTemplateLibrary.makeDAminoAcidTemplates();
		genericTemplateLibrary.loadResEntropy(FileTools.readResource("/config/ResEntropy.dat"));
	}
}
