package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.CapturedIOTest;
import edu.duke.cs.osprey.design.commands.CommandBindingAffinity;
import edu.duke.cs.osprey.design.commands.CommandPartitionFunction;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.*;

class MainTest extends CapturedIOTest {

    @Test
    void noArgumentsPrintsGeneralHelp() {
        String[] argv = {};
        var main = new Main();
        var retVal = main.run(argv);
        assertTrue(mockedOut.toString().contains("Usage"));
        assertEquals(1, retVal);
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

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void commandsRequireADesignFile(String commandName) {
        var noArgsRetval = new Main().run(new String[] {commandName});
        assertNotEquals(Main.Success, noArgsRetval, String.format("Command %s should return an error when no design file specified", commandName));
    }

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void requiresADesignFile(String commandName) {
        String[] argv = {commandName};
        var noArgsRetval = new Main().run(argv);
        assertNotEquals(Main.Success, noArgsRetval);
    }

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void commandFailsWithoutDesignFile(String commandName) {
        String[] argv = { commandName,  "--design" };
        assertNotEquals(Main.Success, new Main().run(argv));
    }

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void designMustBePassedAnExistingFile(String commandName) {
        var pathToDesignFile = "non-existant.yaml";
        String[] argv = { commandName,  "--design", pathToDesignFile };
        assertNotEquals(Main.Success, new Main().run(argv));
        assertTrue(mockedOut.toString().contains("exist"));
    }
}
