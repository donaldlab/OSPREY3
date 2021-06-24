package edu.duke.cs.osprey.design.commands;

import edu.duke.cs.osprey.CapturedIOTest;
import edu.duke.cs.osprey.design.Main;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CommandBindingAffinityTest extends CapturedIOTest {

    private static final String commandName = CommandBindingAffinity.CommandName;

    @Test
    void noArgsPrintsUsage() {
        String[] argv = { commandName };
        new Main().run(argv);
        var out = mockedOut.toString();
        var fmt = "Stdout should contain the word \"%s\"";
        assertTrue(out.contains("Usage"), String.format(fmt, "Usage"));
        assertTrue(out.contains("approximation"), String.format(fmt, "approximation"));
        assertTrue(out.contains("K*"), String.format(fmt, "K*"));
    }
}
