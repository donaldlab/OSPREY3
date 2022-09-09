using Test
using BioSequences
import JSON

include("mut-sig-probs.jl")
using .MutSigProbs

melanoma_probs = """
{
    "GGA" : {
       "GAA" : 0.218572395,
       "GCA" : 0.000958632,
       "GTA" : 0.002185331
    },
    "GAT" : {
       "GCT" : 5.0796e-05,
       "GGT" : 0.000731155,
       "GTT" : 0.001277772
    },
    "GGG" : {
       "GAG" : 0.077262288,
       "GCG" : 0.000476074,
       "GTG" : 0.000736481
    },
    "AAT" : {
       "ACT" : 0.000203184,
       "AGT" : 0.002559217,
       "ATT" : 0.003134217
    },
    "GCA" : {
       "GAA" : 0.001181113,
       "GGA" : 0.00012992,
       "GTA" : 0.001426065
    },
    "GTA" : {
       "GAA" : 9.7361e-05,
       "GCA" : 0.000809435,
       "GGA" : 2.822e-05
    }
}
""" |> JSON.parse

lung_adenocarcinoma_probs = """
{
   "GCT" : {
      "GAT" : 0.010743126,
      "GGT" : 0.002779725,
      "GTT" : 0.005172219
   },
   "GGT" : {
      "GAT" : 0.004323732,
      "GCT" : 0.002596398,
      "GTT" : 0.012852828
   },
   "GTT" : {
      "GAT" : 0.002032707,
      "GCT" : 0.002230074,
      "GGT" : 0.000299724
   }
}
""" |> JSON.parse

function calculate_leu_to_met_mutational_probability()
	# Methionine Codon: ATG
	# 1-step: 
	#	CTG ->(0.010743126) ATG = 0.010743126
	# 2-steps: 
	#	CTG -> (0.002779725) GTG -> (0.004323732) ATG = 1.20187859337e-5
	#	CTG -> (0.005172219) TTG -> (0.002032707) ATG = 1.0513605766832999e-5
	init(lung_adenocarcinoma_probs)
	return calculateFivemerMutationalProbabilities(LongDNASeq("GCTGG"))[AA_M]
end

function calculate_g466e_mutational_probability()
	# G466E with melanoma cancer probabilities
	# Codon TGGAT
	#
	# 1-step: 
	#	GGA ->(0.218572395) GAA = 0.218572395
	# 2-step:
	#	GGA ->(0.000731155) GGG ->(0.077262288) GAG = 5.649070818264e-5
	#	GGA ->(0.218572395) GAA ->(0.002559217) GAG = 0.000559374189014715
	#	GGA ->(0.000958632) GCA ->(0.001181113) GAA = 1.132252717416e-6
	#	GGA ->(0.002185331) GTA ->(9.7361e-05) GAA = 2.12766011491e-7
	init(melanoma_probs)
	return calculateFivemerMutationalProbabilities(LongDNASeq("TGGAT"))[AA_E]
end

@testset "Mutational Signature Probability Calculations" begin
	@test calculate_g466e_mutational_probability() == 0.2191896049159263
	@test calculate_leu_to_met_mutational_probability() == 0.010765658391700534
end
