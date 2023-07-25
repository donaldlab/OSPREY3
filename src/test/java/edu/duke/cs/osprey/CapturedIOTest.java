package edu.duke.cs.osprey;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;

import java.io.*;

public class CapturedIOTest {

    protected OutputStream outputStream, errorStream;
    protected ByteArrayOutputStream mockedOut, mockedError;
    protected InputStream inputStream, mockedIn;

    @BeforeEach
    public void beforeEach() {
        outputStream = System.out;
        inputStream = System.in;
        errorStream = System.err;

        mockedOut = new ByteArrayOutputStream();
        mockedError = new ByteArrayOutputStream();
        mockedIn = new ByteArrayInputStream(new byte[]{});

        System.setOut(new PrintStream(mockedOut));
        System.setIn(mockedIn);
    }

    @AfterEach
    public void afterEach() {
        System.setOut(new PrintStream(outputStream));
        System.setErr(new PrintStream(errorStream));
        System.setIn(inputStream);
    }
}
