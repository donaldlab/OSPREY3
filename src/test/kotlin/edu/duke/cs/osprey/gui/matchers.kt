package edu.duke.cs.osprey.gui

import cuchaz.kludge.tools.abs
import io.kotest.matchers.Matcher
import io.kotest.matchers.MatcherResult


private fun testDoubles(obs: Double?, exp: Double?) =
	if (obs == null || exp == null) {

		MatcherResult(
			obs == exp,
			{ "$obs should be equal to $exp" },
			{ "$obs should not be equal to $exp" }
		)

	} else if (exp.isNaN() || obs.isNaN()) {

		val msg = "NaN values can pass no tests, since by design NaN != NaN"
		MatcherResult(
			false,
			{ msg },
			{ msg }
		)

	} else {
		null
	}

infix fun Double.absolutely(epsilon: Double) = AbsoluteMatcher(this, epsilon)

class AbsoluteMatcher(private val exp: Double?, private val epsilon: Double) : Matcher<Double?> {

	override fun test(value: Double?) =
		testDoubles(value, exp) ?: run {

			// testDoubles garantees these can't be null
			// sadly, the contracts feature isn't smart enough yet to tell this to the compiler
			// so we'll just use runtime assertions instead of compiler guarantees for now
			value!!
			exp!!

			val abserr = (value - exp).abs()
			MatcherResult(
				abserr <= epsilon,
				{ "$value should be equal to $exp, but absolute error $abserr is more than epsilon $epsilon" },
				{ "$value should not be equal to $exp, but absolute error $abserr is within epsilon $epsilon" }
			)
		}
}


infix fun Double.relatively(epsilon: Double) = RelativeMatcher(this, epsilon)

class RelativeMatcher(private val exp: Double?, private val epsilon: Double) : Matcher<Double?> {

	override fun test(value: Double?) =
		testDoubles(value, exp) ?: run {

			// testDoubles garantees these can't be null
			// sadly, the contracts feature isn't smart enough yet to tell this to the compiler
			// so we'll just use runtime assertions instead of compiler guarantees for now
			value!!
			exp!!

			val relerr = (value - exp).abs()/exp.abs()
			MatcherResult(
				relerr <= epsilon,
				{ "$value should be equal to $exp, but relative error $relerr is more than epsilon $epsilon" },
				{ "$value should not be equal to $exp, but relative error $relerr is within epsilon $epsilon" }
			)
		}
}
