
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
    
    public static Options CreateOptions() {

        Options options = new Options();

        Option input = Option.builder("path to input VCF file")
                .numberOfArgs(1)
                .desc("Required, must be uncompressed. Will index it.")
                .argName("input").
                build();

        Option output = Option.builder("path to output VCF file")
                .numberOfArgs(1)
                .desc("Optional, by default prints to STDOUT.")
                .argName("output")
                .build();

        Option DP = Option.builder("integer")
                .numberOfArgs(1)
                .desc("Minimum overall DP of each variant (row) in the input.")
                .argName("DP")
                .build();

        Option sDP = Option.builder("integer")
                .numberOfArgs(1)
                .desc("Minimum DP of each sample.")
                .argName("sDP")
                .build();

        Option missing = Option.builder("double")
                .numberOfArgs(1)
                .desc("Allowed % of missing data in each variant. Example: -missing 5")
                .argName("missing")
                .build();

        Option interval = Option.builder("integer")
                .numberOfArgs(2)
                .desc("Select interval of coordinates, requires -chr . Example: -interval 1000 2000")
                .argName("interval")
                .build();

        Option chr = Option.builder("String")
                .numberOfArgs(1)
                .desc("Select a chromosome of interest. Example: -chr Bd1")
                .argName("chr")
                .build();

        Option pos = Option.builder("integer")
                .numberOfArgs(1)
                .desc("Select a precise position. Example: -pos 1234")
                .argName("pos")
                .build();

        Option sample = Option.builder("String: sample names")
                .desc("Select a subset of samples. Example: -sample sample1,sample2")
                .argName("sample")
                .build();

        Option call = Option.builder("String: sample name")
                .numberOfArgs(1)
                .desc("Requires -pos & -sample. Example: -call sample1 -chr Bd1 -pos 1000")
                .argName("call")
                .build();

        Option bi = new Option("bi", "Select only biallelic genotypes.");

        Option MAF = Option.builder("double")
                .numberOfArgs(1)
                .desc("Minimum allele frequency [0-1]. Example: -MAF 0.05")
                .argName("MAF")
                .build();

        Option maxHet = Option.builder("double")
                .numberOfArgs(1)
                .desc("Maximum % of heterozigous samples. Example: -maxHet 90")
                .argName("maxHet")
                .build();

        Option nr = Option.builder("integer in Kb")
                .numberOfArgs(1)
                .desc("Select top quality variant in selected window size, in Kb. Example: -nr 100")
                .argName("nr")
                .build();

        //Option split = new Option( "split", "Split in two VCF files: i)SNPs and ii)Indels" );
        options.addOption(DP);
        options.addOption(input);
        options.addOption(output);
        options.addOption(sDP);
        options.addOption(missing);
        options.addOption(interval);
        options.addOption(chr);
        options.addOption(pos);
        options.addOption(sample);
        options.addOption(bi);
        options.addOption(MAF);
        options.addOption(maxHet);
        options.addOption(nr);
        options.addOption(call);
        //options.addOption(split);
        return options;
    }

    public static void main(String[] args) throws IOException, ParseException {

        Options options = CreateOptions();

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
    }
}
