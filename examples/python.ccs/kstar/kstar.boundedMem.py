
# First, choose how much memory we want to reserve for the free energy calculation.
# Let's start with 1 GiB for now.
free_energy_mib = 1*1024

# Start osprey with at least that much memory, plus a little more.
import osprey
osprey.start(heapSizeMiB=free_energy_mib + 256)

# Import the compiled conformation spaces module after starting Osprey.
import osprey.ccs


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
    parallelism=osprey.Parallelism(cpuCores=4),

    # The A* node database is what takes up most of the memory in a free energy calculation.
    # So bound the size of that memory using the settings we defined earlier.
    nodeDBMem=free_energy_mib*1024*1024 # convert to bytes
)

# Finally, run K*!
# And wait for it to finish.
# The results will be written to the console.
free_energy_calc.run(kstar)
