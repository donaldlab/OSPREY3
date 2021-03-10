package edu.duke.cs.osprey

import io.kotlintest.*
import io.kotlintest.extensions.RuntimeTagExtension
import io.kotlintest.specs.IntelliMarker
import io.kotlintest.specs.KotlinTestDsl
import kotlin.reflect.jvm.jvmName


/**
 * A KotlinTest spec that allows sharing setup and state between tests
 */
open class SharedSpec(body: SharedSpec.() -> Unit = {}) : IntelliMarker, AbstractSpec() {

	// use a single test instance, so state gets shared rather than re-computed
	override fun isolationMode() = IsolationMode.SingleInstance

	private val testCases = ArrayList<TestCase>()
	override fun testCases(): List<TestCase> = testCases

	private var hasFocus = false

	private fun makeConfig(enabled: Boolean? = null, vararg tags: Tag) =
		TestCaseConfig(
			enabled ?: defaultTestCaseConfig.enabled,
			defaultTestCaseConfig.invocations,
			defaultTestCaseConfig.timeout,
			defaultTestCaseConfig.threads,
			defaultTestCaseConfig.tags + tags,
			defaultTestCaseConfig.extensions
		)

	fun test(name: String, enabled: Boolean? = null, focus: Boolean = false, test: suspend TestContext.() -> Unit) {

		if (focus) {
			hasFocus = true
		}

		testCases.add(TestCase(
			Description.spec(this::class).append(name),
			this,
			test,
			sourceRef(),
			TestType.Test,
			if (focus) {
				makeConfig(enabled, tagFocus)
			} else {
				makeConfig(enabled)
			}
		))
	}

	fun group(name: String, init: Group.() -> Unit) {

		val testCase =
			TestCase(
				Description.spec(this@SharedSpec::class).append(name),
				this@SharedSpec,
				{ /* nothing to do here */ },
				sourceRef(),
				TestType.Container,
				makeConfig(null, tagGroup)
			)
		testCases.add(testCase)

		// run the group code now, at config time
		Group(testCase).init()
	}

	@KotlinTestDsl
	inner class Group(val testCase: TestCase) {

		fun test(name: String, enabled: Boolean? = null, focus: Boolean = false, test: suspend TestContext.() -> Unit) {

			if (focus) {
				hasFocus = true
			}

			testCases.add(TestCase(
				testCase.description.append(name),
				this@SharedSpec,
				test,
				sourceRef(),
				TestType.Test,
				if (focus) {
					makeConfig(enabled, tagFocus)
				} else {
					makeConfig(enabled)
				}
			))
		}

		fun group(name: String, init: Group.() -> Unit) {

			val testCase =
				TestCase(
					testCase.description.append(name),
					this@SharedSpec,
					{ /* nothing to do here */ },
					sourceRef(),
					TestType.Container,
					makeConfig(null, tagGroup)
				)
			testCases.add(testCase)

			// run the group code now, at config time
			Group(testCase).init()
		}
	}

	init {

		// configure (but don't run) all the tests
		body()

		// apply focusing via tags if needed
		if (hasFocus) {
			RuntimeTagExtension.included.apply {
				add(tagFocus)
				add(tagGroup)
			}
		}
	}
}

private object tagFocus: Tag()
private object tagGroup: Tag()

private fun sourceRef(): SourceRef {
	val stack = Throwable().stackTrace
	return stack.dropWhile {
		it.className == SharedSpec::class.jvmName
	}[0].run { SourceRef(lineNumber, fileName) }
}



class TestSharedSpec : SharedSpec({

	var globalCounter = 0

	test("test 1") {
		++globalCounter shouldBe 1
	}

	test("test 2") {
		++globalCounter shouldBe 2
	}

	group("group") {

		var suiteCounter = 0

		test("test 3") {
			++globalCounter shouldBe 3
			++suiteCounter shouldBe 1
		}

		test("test 4") {
			++globalCounter shouldBe 4
			++suiteCounter shouldBe 2
		}
	}
})
