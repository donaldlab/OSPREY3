package osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.json.JSONObject;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Paths;


// NOTE: IntelliJ might show a lot of errors on the imports here
// but IntelliJ is wrong, this code compiles and runs just fine in Gradle
// IntelliJ seems to have quite a few bugs when it comes to gradle buildSrc.

// since IntelliJ can't handle this code, it can't launch the tests from the IDE either
// run the tests with by creating the IntelliJ launcher automatically:
// classpath: osprey.buildsrc.test
// class: osprey.TestJavadocTool
// IntelliJ will complain with red highlights, errors, etc, but it's wrong.
// The test will run. Gradle knows more than IntelliJ does here.

// to recompile the code for tests, you'll have to sync the gradle project yourself
// or find some way to convince IntelliJ to sync gradle as a pre-launch step


public class TestJavadocTool {

	@Test
	public void test()
	throws IOException {

		var out = new JSONObject();

		JavadocTool.runFile(
			"osprey",
			Paths.get("src/test/java/osprey/Test.java").toAbsolutePath(),
			out
		);

		// what's the json look like?
		System.out.println("JSON:\n" + out.toString(4));

		// spot check a few things
		var test = out.getJSONObject("osprey.Test");
		assertThat(test.getJSONObject("type").getString("name"), is("Test"));
		assertThat(test.getJSONObject("type").getString("url"), is("osprey/Test.html"));
		assertThat(test.getJSONObject("javadoc").getString("text"), is("Javadoc comment for Test"));
		var i = field(test, "i");
		assertThat(i.getJSONObject("javadoc").getString("text"), is("javadoc for i"));
		assertThat(i.getString("initializer"), is("1"));
		var stuff = method(test, "stuff(int,float)void");
		assertThat(stuff.getString("name"), is("stuff"));
		assertThat(stuff.getJSONObject("javadoc").getString("text"), is("javadoc for stuff"));
		assertThat(stuff.getString("signature"), is("(int,float)void"));
		assertThat(stuff.getString("returns"), is("void"));
		assertThat(stuff.getJSONArray("args").getJSONObject(0).getString("name"), is("a"));
		assertThat(stuff.getJSONArray("args").getJSONObject(0).getString("type"), is("int"));
		assertThat(stuff.getJSONArray("args").getJSONObject(1).getString("name"), is("b"));
		assertThat(stuff.getJSONArray("args").getJSONObject(1).getString("type"), is("float"));
		var container = field(test, "container");
		assertThat(container.getJSONObject("type").getString("name"), is("Container"));
		assertThat(container.getJSONObject("type").getJSONArray("params").getJSONObject(0).getString("name"), is("java.lang.String"));

		var foo = out.getJSONObject("osprey.Test$Foo");
		assertThat(foo.getJSONObject("type").getString("name"), is("Foo"));
		assertThat(foo.getJSONObject("type").getString("url"), is("osprey/Test.Foo.html"));
		var barIndex = field(foo, "barIndex");
		assertThat(barIndex.getJSONObject("type").getString("name"), is("java.util.Map"));
		assertThat(barIndex.getJSONObject("type").getJSONArray("params").getJSONObject(0).getString("name"), is("java.lang.String"));
		assertThat(barIndex.getJSONObject("type").getJSONArray("params").getJSONObject(1).getString("name"), is("Bar"));
		assertThat(barIndex.getJSONObject("type").getJSONArray("params").getJSONObject(1).getString("url"), is("osprey/Test.Foo.Bar.html"));

		// check some javadoc things too

		// like method arguments
		var javadoc = method(test, "params(int,int)void").getJSONObject("javadoc");
		assertThat(javadoc.getString("text"), is("This method has params"));

		var param = javadoc.getJSONArray("params").getJSONObject(0);
		assertThat(param.getString("text"), is("@param a the a"));
		assertThat(param.getString("name"), is("a"));
		assertThat(param.getString("description"), is("the a"));

		param = javadoc.getJSONArray("params").getJSONObject(1);
		assertThat(param.getString("text"), is("@param b the b"));
		assertThat(param.getString("name"), is("b"));
		assertThat(param.getString("description"), is("the b"));

		// like link parsing
		javadoc = field(test, "j").getJSONObject("javadoc");
		var link = javadoc.getJSONArray("links").getJSONObject(0);
		assertThat(link.getString("text"), is("{@link #container}"));
		assertThat(link.getString("label"), is("container"));
		assertThat(link.getString("signature"), is("container"));
		assertThat(link.getJSONObject("type").getString("name"), is("Test"));

		link = javadoc.getJSONArray("links").getJSONObject(1);
		assertThat(link.getString("text"), is("{@link List}"));
		assertThat(link.getJSONObject("type").getString("name"), is("java.util.List"));

		link = javadoc.getJSONArray("links").getJSONObject(2);
		assertThat(link.getString("text"), is("{@link Foo a foo}"));
		assertThat(link.getJSONObject("type").getString("name"), is("Foo"));
		assertThat(link.getString("label"), is("a foo"));

		link = javadoc.getJSONArray("links").getJSONObject(3);
		assertThat(link.getString("text"), is("{@link #stuff}"));
		assertThat(link.getString("label"), is("stuff"));
		assertThat(link.getString("signature"), is("stuff(int,float)"));

		link = javadoc.getJSONArray("links").getJSONObject(4);
		assertThat(link.getString("text"), is("{@link #stuff(double)}"));
		assertThat(link.getString("label"), is("stuff"));
		assertThat(link.getString("signature"), is("stuff(double)"));

		// and non-standard javadoc extensions
		javadoc = field(test, "k").getJSONObject("javadoc");
		assertThat(javadoc.getJSONArray("citations").getJSONObject(0).getJSONArray("lines").getString(0), is("This is a citation."));
		assertThat(javadoc.getJSONArray("warnings").getJSONObject(0).getString("content"), is("This is a warning."));
		assertThat(javadoc.getJSONArray("notes").getJSONObject(0).getString("content"), is("This is a note."));

		// make sure content gets concatenated correctly
		javadoc = field(test, "m").getJSONObject("javadoc");
		assertThat(javadoc.getString("text"), is("hello {@note hi} world"));

		// check for block tags too
		javadoc = field(test, "blocks").getJSONObject("javadoc");
		assertThat(javadoc.getString("text"), is("This one has block comments."));
		assertThat(javadoc.getJSONArray("notes").getJSONObject(0).getString("content"), is("the note"));

		assertThat(out.has("osprey.Test$Foo$Bar"), is(true));

		// check class-like things too
		var eenumm = out.getJSONObject("osprey.Test$Eenumm");
		assertThat(eenumm.getJSONObject("type").getString("name"), is("Eenumm"));
		assertThat(eenumm.getJSONObject("javadoc").getString("text"), is("This is an eenumm!"));
		assertThat(eenumm.getJSONArray("fields").getJSONObject(0).getString("name"), is("Val1"));
		assertThat(eenumm.getJSONArray("fields").getJSONObject(1).getString("name"), is("Val2"));
	}

	private static JSONObject method(JSONObject c, String id) {
		var i = c.getJSONObject("methodsLut").getInt(id);
		return c.getJSONArray("methods").getJSONObject(i);
	}

	private static JSONObject field(JSONObject c, String id) {
		var i = c.getJSONObject("fieldsLut").getInt(id);
		return c.getJSONArray("fields").getJSONObject(i);
	}
}
