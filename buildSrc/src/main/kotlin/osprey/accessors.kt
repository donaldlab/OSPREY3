package osprey

import org.beryx.runtime.JPackageImageTask
import org.gradle.api.NamedDomainObjectContainer
import org.gradle.api.Project
import org.gradle.api.tasks.SourceSetContainer
import org.gradle.kotlin.dsl.the
import org.jetbrains.kotlin.gradle.dsl.KotlinJvmProjectExtension
import org.beryx.runtime.data.RuntimePluginExtension
import org.gradle.api.Action
import org.gradle.api.NamedDomainObjectProvider
import org.gradle.api.tasks.TaskContainer
import org.gradle.api.tasks.TaskProvider
import org.gradle.kotlin.dsl.named


// gradle doesn't generate some accessors correctly for the Kotlin DSL when used in buildSrc
// so we have to make them ourselves, le sigh =(

// TODO: is there a way to get Gradle to do the right thing here?
//   seems like a complicated issue, see:
//   https://stackoverflow.com/questions/52975515/unresolved-reference-sourcesets-for-gradle-kotlin-dsl
//   https://github.com/gradle/kotlin-dsl-samples/issues/673
//   https://stackoverflow.com/questions/58826817/gradle-pre-compiled-script-plugin-fails-with-expression-cannot-be-invoked-a/58826818

val Project.sourceSets get() =
	the<SourceSetContainer>()

val <T> NamedDomainObjectContainer<T>.main get() =
	getByName("main")

val <T> NamedDomainObjectContainer<T>.test get() =
	getByName("test")


val Project.kotlin: KotlinJvmProjectExtension get() =
	extensions.getByName("kotlin") as KotlinJvmProjectExtension

val Project.runtime: RuntimePluginExtension get() =
    extensions.getByName("runtime") as RuntimePluginExtension

fun Project.runtime(configure: Action<RuntimePluginExtension>): Unit =
	extensions.configure("runtime", configure)

val TaskContainer.jpackageImage: TaskProvider<JPackageImageTask> get() =
	named<JPackageImageTask>("jpackageImage")

operator fun <T> NamedDomainObjectProvider<T>.invoke(action: T.() -> Unit) =
	configure(action)
