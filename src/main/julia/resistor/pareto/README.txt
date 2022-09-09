Pareto Optimization for RESISTOR

This is an implementation of Pareto ranking for the RESISTOR algorithm. It is
written in Julia and developed with Julia 1.7.1. It should probably work with
older Julia versions as well but that has not been tested. If in doubt, use
Julia 1.7.x, where x is any patch release. It was tested on Windows, but
there's nothing operating system-specific about the program and should work on
any Julia-supported environment.

This program expects a CSV and a settings file. The CSV must have row headers,
and the settings file must be in a standard format. You must create a Julia
environment that downloads and contains the requisite packages, and you can
then run the program.

Once the packages for the environment have been downloaded, the program can be
invoked like so:

julia --project=. .\main.jl C:\Users\nsg\cols-to-optimize.csv C:\Users\nsg\pareto-settings.json

The settings file is a JSON file that contains the names of the columns and
whether a column should be maximized or minimized. The content of the settings
file we used for RESISTOR, for example, looked like:

{
	"columns": [
		{
			"name": "Sig Prob",
			"optimize": "maximize"
		},
		{
			"name": "K* ATP",
			"optimize": "maximize"
		},
		{
			"name": "K* Gefitinib",
			"optimize": "minimize"
		},
		{
			"name": "Count",
			"optimize": "maximize"
		}
	]
}

This indicates that there are four columns to optimize over, and that the
numerical values in the column "Sign Prob" should be maximized, "K* ATP"
maximized, "K* Gefitinib", and "Count" maximized.

For reference, the first few rows of the CSV file "cols-to-optimize.csv" look
like this:

Position,WT AA,Mut AA,Sig Prob,K* ATP,K* Gefitinib,Count
718,LEU,PHE,2.47E-04,17.16,-23.15,7
718,LEU,TRP,1.75E-05,17.96,1.76,7
718,LEU,HIP,4.20E-04,18.86,4.5,7
...

This program will write the input CSV to standard out with an addition column,
"Rank", with the Pareto rank of the row given the settings specification.

In writing this program, I was influenced by the implementation of the Pareto
ranking algorithm in the KNIME data analytics software package extension
"Erlwood Knime Open Source Cheminformatics."
