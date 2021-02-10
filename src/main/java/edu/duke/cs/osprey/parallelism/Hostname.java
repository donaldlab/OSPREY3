package edu.duke.cs.osprey.parallelism;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.InetAddress;
import java.net.UnknownHostException;


public class Hostname {

	private static String hostname = null;

	public static String get() {

		if (hostname != null) {
			return hostname;
		}

		// no hostname known yet, try to get it

		// try running a "hostname" command first
		try {
			var process = new ProcessBuilder()
				.command("hostname")
				.start();
			try (var in = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
				hostname = in.readLine();
			}
		} catch (IOException ex) {
			// nope, didn't work
		}

		if (hostname != null) {
			return hostname;
		}

		// finally, try to lookup a hostname for the localhost IP
		try {
			hostname = InetAddress.getLocalHost().getHostName();
		} catch (UnknownHostException ex) {
			// nope, didn't work
		}

		return hostname;
	}

	public static void main(String[] args) {
		System.out.println("hostname: " + get());
	}
}