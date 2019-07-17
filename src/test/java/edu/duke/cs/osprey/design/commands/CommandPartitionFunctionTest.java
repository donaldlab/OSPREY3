package edu.duke.cs.osprey.design.commands;

import edu.duke.cs.osprey.CapturedIOTest;
import edu.duke.cs.osprey.design.Main;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CommandPartitionFunctionTest extends CapturedIOTest {

    private static final String commandName = CommandPartitionFunction.CommandName;

    @Test
    void noArgsPrintsUsage() {
        String[] argv = { commandName };
        new Main().run(argv);
        var out = mockedOut.toString();
        assertTrue(out.contains("Usage"));
        assertTrue(out.contains("partition function"));
    }

    @Test
    void helpCommandExists() {
        String[] argv = { commandName };
        var noArgsRetval = new Main().run(argv);
        var noArgsOut = mockedOut.toString();

        mockedOut.reset();
        var helpArgsRetval = new Main().run(new String[]{commandName, "--help"});
        var helpArgsOut = mockedOut.toString();

        assertNotEquals(noArgsRetval, helpArgsRetval);
        assertEquals(noArgsOut, helpArgsOut);
    }

    @Test
    void requiresADesignFile() {
        String[] argv = { commandName };
        var noArgsRetval = new Main().run(argv);
        assertNotEquals(Main.Success, noArgsRetval);
    }

    @Test
    void commandFailsWithoutDesignFile() {
        String[] argv = { commandName,  "--design" };
        assertNotEquals(Main.Success, new Main().run(argv));
    }

    @Test
    void commandSucceedsWithDesignFile() {
        var pathToDesignFile = getClass().getResource("/test-design.yaml").getPath();
        String[] argv = { commandName,  "--design", pathToDesignFile };
        assertEquals(Main.Success, new Main().run(argv));
    }

    @Test
    void designMustBePassedAnExistingFile() {
        var pathToDesignFile = "non-existant.yaml";
        String[] argv = { commandName,  "--design", pathToDesignFile };
        assertNotEquals(Main.Success, new Main().run(argv));
        assertTrue(mockedOut.toString().contains("exist"));
    }
}
