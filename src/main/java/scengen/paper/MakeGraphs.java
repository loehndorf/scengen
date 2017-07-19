package scengen.paper;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import scengen.distribution.DISTRIBUTION;
import scengen.methods.METHOD;
import scengen.tools.Xorshift;

public class MakeGraphs {

	public static void main(String... args) {
		int numScen = 100;
		Random r = new Xorshift(11);
		List<File> files = new ArrayList<>();
		files.add(makeScen(8, numScen, DISTRIBUTION.Normal, METHOD.QuantizationLearning, r));
		files.add(makeScen(32, numScen, DISTRIBUTION.Normal, METHOD.QuantizationLearning, r));
		files.add(makeScen(8, numScen, DISTRIBUTION.Normal, METHOD.VoronoiCellSampling, r));
		files.add(makeScen(32, numScen, DISTRIBUTION.Normal, METHOD.VoronoiCellSampling, r));
		int dim = 2;
		for (DISTRIBUTION dist : DISTRIBUTION.values())
			for (METHOD method : METHOD.values()) {
				try{
					files.add(makeScen(dim, numScen, dist, method, r));
				}
				catch(IllegalArgumentException e) {
					e.printStackTrace();
				};
			}
		zip(files,"graph_csvs.zip");
		for (File f : files)
			f.delete();
	}
	
	private static File makeScen(int dim, int numScen, DISTRIBUTION dist, METHOD method, Random rand) {
		Map<double[],Double> scen = RunGenerator.makeScen(dim, numScen, dist, method, rand);
		String filename = String.format("scen_%s_%d_%d_%s.csv", dist.name(),dim,numScen,method.name());
		return toCSV(scen, dim, filename);
	}

	private static File toCSV(Map<double[], Double> scen, int dim, String filename) {
		File f = new File(filename);
		FileOutputStream stream;
		BufferedWriter br = null;
		try {
			stream = new FileOutputStream(filename);
			br = new BufferedWriter(new OutputStreamWriter(stream));
			StringBuilder sb = new StringBuilder();
			sb.append("probability");
			for (int i=0; i<dim; i++)
				sb.append(String.format(",value(%d)",i));
			sb.append("\n");
			for (double[] x : scen.keySet()) {
				sb.append(String.format("%f",scen.get(x)));
				for (int i=0; i<dim; i++)
					sb.append(String.format(",%f",x[i]));
				sb.append("\n");
			}
			br.write(sb.toString());
			br.flush();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return f;
	}
	
	public static File zip(List<File> files, String filename) {
	    File zipfile = new File(filename);
	    // Create a buffer for reading the files
	    byte[] buf = new byte[1024];
	    try {
	        // create the ZIP file
	        ZipOutputStream out = new ZipOutputStream(new FileOutputStream(zipfile));
	        // compress the files
	        for(int i=0; i<files.size(); i++) {
	            FileInputStream in = new FileInputStream(files.get(i).getCanonicalFile());
	            // add ZIP entry to output stream
	            out.putNextEntry(new ZipEntry(files.get(i).getName()));
	            // transfer bytes from the file to the ZIP file
	            int len;
	            while((len = in.read(buf)) > 0) {
	                out.write(buf, 0, len);
	            }
	            // complete the entry
	            out.closeEntry();
	            in.close();
	        }
	        // complete the ZIP file
	        out.close();
	        return zipfile;
	    } catch (IOException ex) {
	        System.err.println(ex.getMessage());
	    }
	    return null;
	}

}
