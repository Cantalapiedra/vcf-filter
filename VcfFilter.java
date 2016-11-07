
/**
 *
 * @author Eduardo Candeal 2016
 *
 * Refactor CPCantalapiedra 2016
 *
 */
import java.io.IOException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.MissingOptionException;

public class VcfFilter {

    private static final boolean debuggin = true;
    
    /**
     * 
     * Method which sets up the Options available
     * from the command line
     * 
     */
    public static Options CreateOptions() {

        Options options = new Options();
        
        // Input file
        Option input = Option.builder("i")
                .longOpt("input")
                .desc("VCF/BCF type input file.")
                .argName("FILE PATH")
                .hasArg(true)
                .numberOfArgs(1)
                .type(String.class)
                .required(true)
                .build();
        options.addOption(input);
        
        // Input file type
        Option inputType = Option.builder("b")
                .longOpt("bcf")
                .desc("Input file format is BCF (default input file format is VCF).")
                .hasArg(false)
                .build();
        options.addOption(inputType);
        
        // Output file
        Option output = Option.builder("o")
                .longOpt("output")
                .desc("Output file (stdout if absent).")
                .argName("FILE PATH")
                .numberOfArgs(1)
                .type(String.class)
                .build();
        options.addOption(output);
        
        // Subsetting (filtering by position) options
        //
        // Chromosome
        Option chr = Option.builder("chr")
                .longOpt("chromosome")
                .desc("Chromosome to be output. Example: -chr Bd1.")
                .hasArg(true)
                .numberOfArgs(1)
                .type(Number.class)
                .argName("Chr name")
                .build();
        options.addOption(chr);
        
        // Single position to be output
        Option pos = Option.builder("pos")
                .longOpt("position")
                .desc("Select a precise position. Requires -chr. Example: -pos 1234")
                .hasArg(true)
                .numberOfArgs(1)
                .type(Number.class)
                .argName("int")
                .build();
        options.addOption(pos);
        
        // Interval of positions, or range.
        Option interval = Option.builder("r")
                .longOpt("range")
                .desc("Interval of positions to be output. Requires -chr. Example: -interval 1000 2000")
                .hasArg(true)
                .numberOfArgs(2)
                .type(Number.class)
                .argName("int int")
                .build();
        options.addOption(interval);
        
        // General depth of variant
        Option DP = Option.builder("dp")
                .longOpt("depth")
                .desc("Minimum overall DP of each variant (row) in the input.")
                .hasArg(true)
                .numberOfArgs(1)
                .type(Number.class)
                .argName("int")
                .build();
        options.addOption(DP);
        
        // Select specific samples to be output, discarding the others
        //
        // By sample name
        Option sample = Option.builder("s")
                .longOpt("samples")
                .desc("Select a subset of samples. Example: -sample sample1,sample2")
                .hasArgs()
                .numberOfArgs(Option.UNLIMITED_VALUES)
                .valueSeparator(',')
                .type(String.class)
                .argName("sample1,sample2,...,sampleN")
                .build();
        options.addOption(sample);
        
        // Missing samples
        Option missing = Option.builder("m")
                .longOpt("missing")
                .desc("Allowed % of missing data in each variant. Example: -missing 5")
                .hasArg(true)
                .numberOfArgs(1)
                .type(Number.class)
                .argName("int")
                .build();
        options.addOption(missing);
        
        // Minimum allele frequency
        Option MAF = Option.builder("MAF")
                .longOpt("min_allele_freq")
                .desc("Minimum allele frequency [0-1]. Example: -MAF 0.05")
                .hasArg(true)
                .numberOfArgs(1)
                .type(Number.class)
                .argName("decimal")
                .build();
        options.addOption(MAF);
        
        // Max heterozygots
        Option maxHet = Option.builder("h")
                .longOpt("max_het")
                .desc("Maximum % of heterozigous samples. Example: -maxHet 90")
                .hasArg(true)
                .numberOfArgs(1)
                .type(Number.class)
                .argName("decimal")
                .build();
        options.addOption(maxHet);
        
        // Minimum depth of each sample
        Option sDP = Option.builder("sdp")
                .longOpt("sample_depth")
                .desc("Minimum DP of each sample.")
                .hasArg(true)
                .numberOfArgs(1)
                .type(Number.class)
                .argName("int")
                .build();
        options.addOption(sDP);
        
        Option call = Option.builder("call")
                .numberOfArgs(1)
                .desc("Requires -pos & -sample. Example: -call sample1 -chr Bd1 -pos 1000")
                .argName("String: sample name")
                .build();
        options.addOption(call);

        Option bi = new Option("bi", "Select only biallelic genotypes.");
        options.addOption(bi);        
        
        Option nr = Option.builder("nr")
                .numberOfArgs(1)
                .desc("Select top quality variant in selected window size, in Kb. Example: -nr 100")
                .argName("integer in Kb")
                .build();
        options.addOption(nr);
        
        //Option split = new Option( "split", "Split in two VCF files: i)SNPs and ii)Indels" );
        //options.addOption(split);
        
        return options;
    }

    public static void main(String[] args) throws Exception {
        
        Options options = CreateOptions();
        
        try {
            CommandLine cmd = new DefaultParser().parse(options, args);

            if (args.length < 1) {
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp(" ", options);

                System.out.println("\nThis utility builds on HTSJDK and handles "
                        + "VCF versions supported there, currently v4.2.");
                System.out.println("Eduardo Candeal, with help from "
                        + "Carlos P Cantalapiedra and Bruno Contreras-Moreira\nEEAD-CSIC 2016");
            } else {
                VcfController vcfctrl = new VcfController(cmd);
                vcfctrl.run();
            }
        
        } catch (ParseException pe) {
            System.out.println("Error when parsing program arguments:");
            System.out.println(pe.getMessage());
            System.out.println();
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", options);
            
            if (debuggin) throw pe;
            
        } catch (IOException ie) {
            System.out.println("I/O error:");
            System.out.println(ie.getMessage());
            if (debuggin) throw ie;
            
        } catch (IllegalArgumentException iae) {
            System.out.println("Error when reading program arguments:");
            System.out.println(iae.getMessage());
            System.out.println();
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", options);
            
            if (debuggin) throw iae;
            
        } catch (Exception e) {
            System.out.println("Error:");
            System.out.println(e.getMessage());
            if (debuggin) throw e;
            
        } finally {
            System.out.println();
        }
    }
}
