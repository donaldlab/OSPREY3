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

		var code = JavadocTool.runFile(Paths.get("src/test/java/build/Test.java"), out);
		assertThat(code, is(0));

		// spot check a few things
		var test = out.getJSONObject("build.Test");
		assertThat(test.getString("name"), is("Test"));
		assertThat(test.getString("javadoc"), is("Javadoc comment for Test"));
		var i = test.getJSONObject("fields").getJSONObject("i");
		assertThat(i.getString("javadoc"), is("javadoc for i"));
		assertThat(i.getString("initializer"), is("1"));
		var stuff = test.getJSONObject("methods").getJSONObject("stuff");
		assertThat(stuff.getString("javadoc"), is("javadoc for stuff"));
		assertThat(stuff.getString("signature"), is("(int,float)void"));
		assertThat(stuff.getString("returns"), is("void"));

		assertThat(out.has("build.Test$Foo"), is(true));
		assertThat(out.has("build.Test$Foo$Bar"), is(true));
	}
}
