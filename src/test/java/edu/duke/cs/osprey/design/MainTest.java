package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.CapturedIOTest;
import edu.duke.cs.osprey.design.commands.CommandBindingAffinity;
import edu.duke.cs.osprey.design.commands.CommandPartitionFunction;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.*;

class MainTest extends CapturedIOTest {

    private final String ke07ApoPath = getClass().getResource("/ke07-apo.yaml").getPath();
    private final String testDesignPath = getClass().getResource("/test-design.yaml").getPath();

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
    void commandSucceedsWithDesignFile(String commandName) {
        String[] argv = { commandName,  "--design", testDesignPath };
        assertEquals(Main.Success, new Main().run(argv));
    }

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void designMustBePassedAnExistingFile(String commandName) {
        var pathToDesignFile = "non-existant.yaml";
        String[] argv = { commandName,  "--design", pathToDesignFile };
        assertNotEquals(Main.Success, new Main().run(argv));
        assertTrue(mockedOut.toString().contains("exist"));
    }

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void epsilonShouldBeAnOption(String commandName) {
        String[] argv = { commandName, "--design", ke07ApoPath, "--epsilon", "0.7" };
        assertEquals(Main.Success, new Main().run(argv));
    }

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void epsilonRequiresAnArgument(String commandName) {
        String[] argv = { commandName, "--design", ke07ApoPath, "--epsilon" };
        assertNotEquals(Main.Success, new Main().run(argv));
    }

    @ParameterizedTest
    @ValueSource(strings = {CommandBindingAffinity.CommandName, CommandPartitionFunction.CommandName})
    void epsilonMustBeAFloat(String commandName) {
        String[] argv = { commandName, "--design", ke07ApoPath, "--epsilon", "one-point-oh" };
        assertNotEquals(Main.Success, new Main().run(argv));
    }
}
