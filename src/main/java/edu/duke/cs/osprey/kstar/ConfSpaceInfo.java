package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.AutoCloseableNoEx;

import java.io.File;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;

public class ConfSpaceInfo {

    private final KStarSettings settings;
    public final ConfSpaceIteration confSpace;
    public final ConfSpaceType type;
    public ConfEnergyCalculator confEcalc = null;
    public final String id;

    public final Map<Sequence, PartitionFunction.Result> pfuncResults = new HashMap<>();

    public File confDBFile = null;
    public PfuncFactory pfuncFactory = null;

    public interface PfuncFactory {
        PartitionFunction make(RCs rcs);
    }
    private ConfDB confDB = null;

    public ConfSpaceInfo(KStarSettings settings, ConfSpaceIteration confSpace, ConfSpaceType type) {
        this.settings = settings;
        this.confSpace = confSpace;
        this.type = type;
        this.id = type.name().toLowerCase();

        confDBFile = new File(String.format(settings.confDBPattern, id));
    }

    public void check() {
        if (confEcalc == null) {
            throw new InitException(type, "confEcalc");
        }
        if (pfuncFactory == null) {
            throw new InitException(type, "pfuncFactory");
        }
    }

    public AutoCloseableNoEx openConfDB() {
        if (confDBFile != null) {
            /*
            if (!settings.resume) {
                confDBFile.delete();
            }
             */
            confDB = new ConfDB(confSpace, confDBFile);
        }
        return () -> {
            if (confDB != null) {
                confDB.close();
                confDB = null;
            }
        };
    }

    public void clear() {
        pfuncResults.clear();
    }

    public PartitionFunction.Result calcPfunc(TaskExecutor.ContextGroup ctxGroup, Sequence globalSequence, BigDecimal stabilityThreshold) {

        Sequence sequence = globalSequence.filter(confSpace.seqSpace());

        // check the cache first
        PartitionFunction.Result result = pfuncResults.get(sequence);
        if (result != null) {
            return result;
        }

        // cache miss, need to compute the partition function

        // compute the partition function
        PartitionFunction pfunc = makePfunc(ctxGroup, sequence);
        pfunc.setStabilityThreshold(stabilityThreshold);
        if (settings.pfuncTimeout != null) {
            pfunc.compute(settings.pfuncTimeout);
        } else {
            pfunc.compute(settings.maxNumConfs);
        }

        // save the result
        result = pfunc.makeResult();
        pfuncResults.put(sequence, result);

        /* HACKHACK: we're done using the A* tree, pfunc, etc
            and normally the garbage collector will clean them up,
            along with their off-heap resources (e.g. TPIE data structures).
            Except the garbage collector might not do it right away.
            If we try to allocate more off-heap resources before these get cleaned up,
            we might run out. So poke the garbage collector now and try to get
            it to clean up the off-heap resources right away.
        */
        Runtime.getRuntime().gc();

        // newer JVMs have more concurrent garbage collectors
        // give it a little time to finish cleaning up the pfunc
        try {
            Thread.sleep(10);
        } catch (InterruptedException ex) {
            throw new RuntimeException(ex);
        }

        return result;
    }

    public PartitionFunction makePfunc(TaskExecutor.ContextGroup ctxGroup, Sequence seq) {

        RCs rcs = seq.makeRCs(confSpace);

        PartitionFunction pfunc = pfuncFactory.make(rcs);

        pfunc.setReportProgress(settings.showPfuncProgress);
        if (settings.useExternalMemory) {
            PartitionFunction.WithExternalMemory.setOrThrow(pfunc, true, rcs);
        }
        if (confDB != null) {
            PartitionFunction.WithConfDB.cast(pfunc).setConfDB(confDB, seq);
        }

        pfunc.setInstanceId(type.ordinal());
        pfunc.putTaskContexts(ctxGroup);

        pfunc.init(settings.epsilon);

        return pfunc;
    }
}
