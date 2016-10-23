package edu.duke.cs.osprey.tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

public class ResourceExtractor {
	
	public static File extract(String resourcePath)
	throws IOException {
		return extract(resourcePath, ResourceExtractor.class);
	}
	
	public static File extract(String resourcePath, Class<?> relativeTo)
	throws IOException {
		
		try (InputStream in = relativeTo.getResourceAsStream(resourcePath)) {
			
			if (in == null) {
				throw new FileNotFoundException("can't find native library at: " + resourcePath);
			}
			
			// copy the library to a temp file
			File file = File.createTempFile(FilenameUtils.getBaseName(resourcePath) + ".", "." + FilenameUtils.getExtension(resourcePath));
			file.deleteOnExit();
			try (OutputStream out = new FileOutputStream(file)) {
				IOUtils.copy(in, out);
			}
			
			return file;
		}
	}
}
