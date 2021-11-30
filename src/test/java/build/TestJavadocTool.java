package build;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.json.JSONObject;
import org.junit.Test;

import java.nio.file.Paths;


public class TestJavadocTool {

	@Test
	public void test() {

		var out = new JSONObject();

		var code = JavadocTool.runFile(
			"build",
			Paths.get("src/test/java/build/Test.java"),
			out
		);
		assertThat(code, is(0));

		// what's the json look like?
		System.out.println("JSON:\n" + out.toString(4));

		// spot check a few things
		var test = out.getJSONObject("build.Test");
		assertThat(test.getJSONObject("type").getString("name"), is("Test"));
		assertThat(test.getJSONObject("type").getString("url"), is("build/Test.html"));
		assertThat(test.getJSONObject("javadoc").getString("text"), is("Javadoc comment for Test"));
		var i = test.getJSONObject("fields").getJSONObject("i");
		assertThat(i.getJSONObject("javadoc").getString("text"), is("javadoc for i"));
		assertThat(i.getString("initializer"), is("1"));
		var stuff = test.getJSONObject("methods").getJSONObject("stuff(int,float)void");
		assertThat(stuff.getString("name"), is("stuff"));
		assertThat(stuff.getJSONObject("javadoc").getString("text"), is("javadoc for stuff"));
		assertThat(stuff.getString("signature"), is("(int,float)void"));
		assertThat(stuff.getString("returns"), is("void"));
		assertThat(stuff.getJSONArray("args").getJSONObject(0).getString("name"), is("a"));
		assertThat(stuff.getJSONArray("args").getJSONObject(0).getString("type"), is("int"));
		assertThat(stuff.getJSONArray("args").getJSONObject(1).getString("name"), is("b"));
		assertThat(stuff.getJSONArray("args").getJSONObject(1).getString("type"), is("float"));
		var container = test.getJSONObject("fields").getJSONObject("container");
		assertThat(container.getJSONObject("type").getString("name"), is("Container"));
		assertThat(container.getJSONObject("type").getJSONArray("params").getJSONObject(0).getString("name"), is("java.lang.String"));

		var foo = out.getJSONObject("build.Test$Foo");
		assertThat(foo.getJSONObject("type").getString("name"), is("Foo"));
		assertThat(foo.getJSONObject("type").getString("url"), is("build/Test.Foo.html"));
		var barIndex = foo.getJSONObject("fields").getJSONObject("barIndex");
		assertThat(barIndex.getJSONObject("type").getString("name"), is("java.util.Map"));
		assertThat(barIndex.getJSONObject("type").getJSONArray("params").getJSONObject(0).getString("name"), is("java.lang.String"));
		assertThat(barIndex.getJSONObject("type").getJSONArray("params").getJSONObject(1).getString("name"), is("Bar"));
		assertThat(barIndex.getJSONObject("type").getJSONArray("params").getJSONObject(1).getString("url"), is("build/Test.Foo.Bar.html"));

		// check some javadoc things too

		// like method arguments
		var javadoc = test.getJSONObject("methods").getJSONObject("params(int,int)void").getJSONObject("javadoc");
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
		javadoc = test.getJSONObject("fields").getJSONObject("j").getJSONObject("javadoc");
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
		javadoc = test.getJSONObject("fields").getJSONObject("k").getJSONObject("javadoc");
		assertThat(javadoc.getJSONArray("citations").getJSONObject(0).getJSONArray("lines").getString(0), is("This is a citation."));
		assertThat(javadoc.getJSONArray("warnings").getJSONObject(0).getString("content"), is("This is a warning."));
		assertThat(javadoc.getJSONArray("notes").getJSONObject(0).getString("content"), is("This is a note."));

		// make sure content gets concatenated correctly
		javadoc = test.getJSONObject("fields").getJSONObject("m").getJSONObject("javadoc");
		assertThat(javadoc.getString("text"), is("hello {@note hi} world"));

		// check for block tags too
		javadoc = test.getJSONObject("fields").getJSONObject("blocks").getJSONObject("javadoc");
		assertThat(javadoc.getString("text"), is("This one has block comments."));
		assertThat(javadoc.getJSONArray("notes").getJSONObject(0).getString("content"), is("the note"));

		assertThat(out.has("build.Test$Foo$Bar"), is(true));
	}
}
