import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertNotNull;

public class ForceFieldParamsTest {

    @Test
    void canLoadForcefieldParams() {
        var params = new ForcefieldParams();
        assertNotNull(params.solvationForcefield);
    }
}
