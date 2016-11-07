
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

public class VcfFilter {

    private static boolean debuggin = false;
    
    public static Options CreateOptions() {

        Options options = new Options();
        
        Option input = Option.builder("input")
                .numberOfArgs(1)
                .desc("Required, must be uncompressed. Will index it.")
                .argName("path to input VCF file").
                build();
        options.addOption(input);

        Option output = Option.builder("output")
                .numberOfArgs(1)
                .desc("Optional, by default prints to STDOUT.")
                .argName("path to output VCF file")
                .build();
        options.addOption(output);

        Option DP = Option.builder("DP")
                .numberOfArgs(1)
                .desc("Minimum overall DP of each variant (row) in the input.")
                .argName("integer")
                .build();
        options.addOption(DP);

        Option sDP = Option.builder("sDP")
                .numberOfArgs(1)
                .desc("Minimum DP of each sample.")
                .argName("integer")
                .build();
        options.addOption(sDP);

        Option missing = Option.builder("missing")
                .numberOfArgs(1)
                .desc("Allowed % of missing data in each variant. Example: -missing 5")
                .argName("double")
                .build();
        options.addOption(missing);

        Option interval = Option.builder("interval")
                .numberOfArgs(2)
                .desc("Select interval of coordinates, requires -chr . Example: -interval 1000 2000")
                .argName("integer")
                .build();
        options.addOption(interval);

        Option chr = Option.builder("chr")
                .numberOfArgs(1)
                .desc("Select a chromosome of interest. Example: -chr Bd1")
                .argName("String")
                .build();
        options.addOption(chr);

        Option pos = Option.builder("pos")
                .numberOfArgs(1)
                .desc("Select a precise position. Example: -pos 1234")
                .argName("integer")
                .build();
        options.addOption(pos);

        Option sample = Option.builder("sample")
                .desc("Select a subset of samples. Example: -sample sample1,sample2")
                .argName("String: sample names")
                .build();
        options.addOption(sample);

        Option call = Option.builder("call")
                .numberOfArgs(1)
                .desc("Requires -pos & -sample. Example: -call sample1 -chr Bd1 -pos 1000")
                .argName("String: sample name")
                .build();
        options.addOption(call);

        Option bi = new Option("bi", "Select only biallelic genotypes.");
        options.addOption(bi);        
        
        Option MAF = Option.builder("MAF")
                .numberOfArgs(1)
                .desc("Minimum allele frequency [0-1]. Example: -MAF 0.05")
                .argName("double")
                .build();
        options.addOption(MAF);

        Option maxHet = Option.builder("maxHet")
                .numberOfArgs(1)
                .desc("Maximum % of heterozigous samples. Example: -maxHet 90")
                .argName("double")
                .build();
        options.addOption(maxHet);

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

                System.out.println("\nThis utility builds on HTSJDK and handles VCF versions supported there, currently v4.2.");
                System.out.println("Eduardo Candeal, with help from Carlos P Cantalapiedra and Bruno Contreras-Moreira\nEEAD-CSIC 2016");
            } else {
                VcfController vcfctrl = new VcfController();
                vcfctrl.ReadOptions(cmd);
            }
        
        } catch (ParseException pe) {
            System.out.println("Error when parsing program arguments:");
            System.out.println(pe.getMessage());
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
