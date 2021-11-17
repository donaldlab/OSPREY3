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
		assertThat(test.getString("javadoc"), is("Javadoc comment for Test"));
		var i = test.getJSONObject("fields").getJSONObject("i");
		assertThat(i.getString("javadoc"), is("javadoc for i"));
		assertThat(i.getString("initializer"), is("1"));
		var stuff = test.getJSONObject("methods").getJSONObject("stuff(int,float)void");
		assertThat(stuff.getString("name"), is("stuff"));
		assertThat(stuff.getString("javadoc"), is("javadoc for stuff"));
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

		assertThat(out.has("build.Test$Foo$Bar"), is(true));
	}
}
