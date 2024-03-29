package com.github.lukfor;

import java.util.concurrent.Callable;

import com.github.lukfor.commands.ReportCommand;
import com.github.lukfor.commands.VersionCommand;

import picocli.CommandLine;
import picocli.CommandLine.Command;

public class App {

	public static final String APP = "gwas-report";

	public static final String VERSION = "0.3.1";

	public static final String URL = "https://github.com/lukfor/gwas-report";

	public static final String COPYRIGHT = "(c) 2023 Lukas Forer";

	public static String[] ARGS = new String[0];

	public static void main(String[] args) {

		System.err.println();
		System.err.println(APP + " " + VERSION);
		if (URL != null && !URL.isEmpty()) {
			System.err.println(URL);
		}
		if (COPYRIGHT != null && !COPYRIGHT.isEmpty()) {
			System.err.println(COPYRIGHT);
		}
		System.err.println();

		ARGS = args;

		int exitCode = new CommandLine(new DefaultCommand()).execute(args);
		System.exit(exitCode);

	}

	@Command(name = App.APP, version = App.VERSION, subcommands = { ReportCommand.class, VersionCommand.class })
	public static class DefaultCommand implements Callable<Integer> {

		@Override
		public Integer call() throws Exception {
			// TODO Auto-generated method stub
			return null;
		}

	}

}
