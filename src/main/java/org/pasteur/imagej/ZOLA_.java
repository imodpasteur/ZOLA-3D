

package org.pasteur.imagej;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * A template for processing each pixel of either
 * GRAY8, GRAY16, GRAY32 or COLOR_RGB images.
 *
 * @author Benoit Lelandais
 */
public class ZOLA_ implements PlugInFilter {
	

	@Override
	public int setup(String arg, ImagePlus imp) {
		
		return 0;
	}

	@Override
	public void run(ImageProcessor ip) {
		
	}

	

	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads
	 * an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = ZOLA_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		//String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length()- "/classes".length());
                System.setProperty("plugins.dir", pluginsDir);
		// start ImageJ
		new ImageJ();

		
	}
}
