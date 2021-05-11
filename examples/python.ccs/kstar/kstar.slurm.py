
# First, choose how much memory we want to reserve for the free energy calculation.
# Let's start with 1 GiB for now, but you may want a (much) larger amount for your designs.
free_energy_mib = 1*1024

# Then, pick now much memory we want to allocate for the rest of JVM heap.
# If you're seeing errors like this: "java.lang.OutOfMemoryError: Java heap space",
# try increasing the JVM heap size.
heap_mib = 256

# Then, pick how much memory we want to allocate for the rest of the whole process.
# If SLURM is canceling your jobs for using too much memory, try increasing this.
process_mib = 256

# Start osprey
import osprey
osprey.start(heapSizeMiB=free_energy_mib + heap_mib)

# Import the compiled conformation spaces module after starting Osprey.
import osprey.ccs, osprey.slurm

# Then, pick the cluster, CPU, and GPU resources we want to allocate.
num_nodes = 2
parallelism = osprey.Parallelism(cpuCores=2)

# Launch the job on the SLURM cluster using our resource requirements.
# Make sure to run this script on the login node of your SLURM cluster!
if not osprey.slurm.is_slurm:
    osprey.slurm.launch(num_nodes, parallelism, free_energy_mib + heap_mib + process_mib)

    # after launching the SLURM job, stop processing this script on the login node
    exit(0)

# If we made it this far, we're running on a SLURM compute node and ready to do the processing! =D


kstar = osprey.ccs.kstarBoundedMem(

    # Load the conformation spaces for your design.
    # These conformation spaces were prepared in a previous step.
    # See the 'conf space prep' example.
    complex=osprey.ccs.loadConfSpace("complex.ccsx"),
    design=osprey.ccs.loadConfSpace("ptpase.ccsx"),
    target=osprey.ccs.loadConfSpace("hepes.ccsx"),

    # Search for up to double mutants.
    # Or set to None to look at all combinations of mutations.
    maxSimultaneousMutations=2,

    # Set the precision for calculated free energy values.
    # ie, the width of the interval of lower and upper bounds on the free energy.
    # A smaller width is more precise, but can take longer to calculate.
    # A larger width is less precise, but can take less time to calculate.
    gWidthMax=0.1,

    # Designs that stabilize the bound state but also destabilize the unbound state may not be very useful.
    # Filter out sequences whose unbound states have been destabilized by setting a stability threshold.
    # The filter works by ignoring sequences whose unbound state free energies
    # are too far above the wild type free energy.
    stabilityThreshold=5.0, # in kcal/mol

    # Periodically write out a PDB ensemble of the 5 lowest energy structures for a sequence to the current folder.
    ensembleTracking=[5, '.'],

    # Show progress for each free energy calculation for each sequence for each state.
    reportStateProgress=True
)


# Configure a free energy calculator that uses bounded memory.
free_energy_calc = osprey.ccs.freeEnergyCalc(

    # Reference the multi-state conformation space we defined for the K* design.
    confSpace=kstar.confSpace,

    # More cores is more faster!
    parallelism=parallelism,

    # The A* node database is what takes up most of the memory in a free energy calculation.
    # So bound the size of that memory using the settings we defined earlier.
    nodeDBMem=free_energy_mib*1024*1024 # convert to bytes
)

# Finally, run K*!
# And wait for it to finish.
# The results will be written to the console.
free_energy_calc.run(kstar)
