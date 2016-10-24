package edu.duke.cs.osprey.tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

public class ResourceExtractor {
	
	public static File extract(URL url)
	throws IOException {
		
		if (url == null) {
			throw new FileNotFoundException("url can't be null");
		}
		
		try (InputStream in = url.openStream()) {
			
			// copy the library to a temp file
			String filename = url.getFile();
			File file = File.createTempFile(FilenameUtils.getBaseName(filename) + ".", "." + FilenameUtils.getExtension(filename));
			file.deleteOnExit();
			try (OutputStream out = new FileOutputStream(file)) {
				IOUtils.copy(in, out);
			}
			
			return file;
		}
	}
}
