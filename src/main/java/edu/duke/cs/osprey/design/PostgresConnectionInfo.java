package edu.duke.cs.osprey.design;

public class PostgresConnectionInfo {
    public final String username;
    public final String password;
    public final String connectionString;

    public PostgresConnectionInfo(String username, String password, String connectionString) {
        this.username = username;
        this.password = password;
        this.connectionString = connectionString;
    }
}
