package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.CapturedIOTest;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class MainTest extends CapturedIOTest {

    @Test
    void noArgumentsPrintsGeneralHelp() {
        String[] argv = {};
        var main = new Main();
        var retVal = main.run(argv);
        assertTrue(mockedOut.toString().contains("Usage"));
        assertEquals(0, retVal);
    }

    @Test
    void askingForHelpPrintsHelp() {
        String[] argv = {"--help"};
        new Main().run(argv);
        assertTrue(mockedOut.toString().contains("Usage"));
    }

    @Test
    void aCommandThatDoesNotExistReturnsNonZero() {
        String[] argv = {"not-real"};
        assertNotEquals(0, new Main().run(argv));
    }

    @Test
    void stabilityIsACommand() {
        String[] argv = {"stability", "--help"};
        assertEquals(0, new Main().run(argv));
    }

    @Test
    void affinityIsACommand() {
        String[] argv = {"affinity", "--help"};
        assertEquals(0, new Main().run(argv));
    }
}
