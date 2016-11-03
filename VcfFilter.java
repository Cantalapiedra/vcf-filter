
/**
 *
 * @author Eduardo Candeal 2016
 */

import htsjdk.tribble.index.AbstractIndex;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.io.IOException;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFCodec;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Set;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;

public class VcfFilter {

    public int num;
    public String ifile;
    public String pathin = "";
    public String pathout = "";

    public VCFFileReader VCFreader;

    public static CommandLine cmd;

    public Split split;
    public int DPG;
    public int DPS;
    public String variante;
    public String sample;
    public String crom;
    public int position;
    public int posfirstint;
    public int possecondint;
    public double minHet;
    public int minGQ;
    public int Size;
    public Set set;

    public int total;
    public double frecrar;
    public double NData;
    public VariantContextWriter vcfwriter;
    public EnumSet<Options> DEFAULT_OPTIONS
            = EnumSet.of(Options.INDEX_ON_THE_FLY,
                    Options.ALLOW_MISSING_FIELDS_IN_HEADER,
                    Options.WRITE_FULL_FORMAT_FIELD);

    public void CreateVCF() throws IOException {
        VCFreader = new VCFFileReader(new File(pathin));
        if (pathout.isEmpty()) {
            vcfwriter = VariantContextWriterFactory.createVcf(null, System.out,
                    VCFreader.getFileHeader().getSequenceDictionary(), DEFAULT_OPTIONS);
        } else {
            FileOutputStream outputstream = new FileOutputStream(new File(pathout));
            vcfwriter = VariantContextWriterFactory.createVcf(new File(pathout), outputstream,
                    VCFreader.getFileHeader().getSequenceDictionary(), DEFAULT_OPTIONS);

        }

        vcfwriter.writeHeader(VCFreader.getFileHeader());
    }

    public static void writeTribbleIndex(Index idx, String idxFile) throws IOException {
        LittleEndianOutputStream stream = null;
        try {
            stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
            idx.write(stream);
        } catch (Exception e) {

            // Delete output file as its probably corrupt
            File tmp = new File(idxFile);
            if (tmp.exists()) {
                tmp.delete();
            }
            throw e;
        } finally {
            if (stream != null) {
                stream.close();
            }
        }
    }

    public void createidx() throws IOException {
        int binSize = 16000;
        File inputFile = new File(pathin);
        VCFCodec codec = new VCFCodec();

        String idxFile = pathin + ".idx";

        AbstractIndex idx = IndexFactory.createLinearIndex(inputFile, codec, binSize);

        //try{
        writeTribbleIndex(idx, idxFile);
        //} catch (Exception e){}
    }

    public void ReadOptions() throws IOException {

        String variante = "";

        if (cmd.hasOption("input")) {
            pathin = cmd.getOptionValues("input")[0];
            createidx();
        }

        if (pathin.isEmpty()) {
            /* VCFreader = new VCFFileReader(System.in); // does not work as 
            VCFFileReader constructor expects a java.io.File file, not a stream */
        } else {
            VCFreader = new VCFFileReader(new File(pathin));
        }

        if (cmd.hasOption("output")) {
            pathout = cmd.getOptionValues("output")[0];
        }

        if (cmd.hasOption("DP")) {
            DPG = Integer.parseInt(cmd.getOptionValues("DP")[0]);
            CreateVCF();
            VcfProcedures.MinDPGen(VCFreader, DPG, this.vcfwriter, variante);
        }
        if (cmd.hasOption("sDP")) {
            DPS = Integer.parseInt(cmd.getOptionValues("sDP")[0]);
            CreateVCF();
            VcfProcedures.MinDPSample(VCFreader, DPS, this.vcfwriter, variante);
        }

        if (cmd.hasOption("missing")) {
            NData = Integer.parseInt(cmd.getOptionValues("missing")[0]);
            CreateVCF();
            VcfProcedures.MissingData(VCFreader, num, this.vcfwriter, variante);
        }

        if (cmd.hasOption("chr")) {
            crom = cmd.getOptionValues("chr")[0];
        }

        if (cmd.hasOption("interval")) {
            posfirstint = Integer.parseInt(cmd.getOptionValues("interval")[0]);
            possecondint = Integer.parseInt(cmd.getOptionValues("interval")[1]);
            CreateVCF();
            VcfProcedures.VariantInter(VCFreader, posfirstint, possecondint,
                    this.crom, vcfwriter, variante);
        }

        if (cmd.hasOption("pos")) {
            position = Integer.parseInt(cmd.getOptionValues("pos")[0]);
        }

        if (cmd.hasOption("sample")) {

            this.set = new HashSet();
            String[] samples2 = cmd.getOptionValues("sample")[0].split(",");
            System.out.println(samples2.length);
            for (int i = 0; i < samples2.length; i++) {
                System.out.println(samples2[i]);
                this.set.add(samples2[i]);
            }
            VcfProcedures.SelectGenotype(VCFreader, this.set, this.pathout,
                    this.vcfwriter, this.DEFAULT_OPTIONS);

        }

        if (cmd.hasOption("call")) {

            sample = cmd.getOptionValues("call")[0];
            VcfProcedures.FindSamVar(VCFreader, position, num,
                    sample, variante);
            System.out.println(variante);

        }
        if (cmd.hasOption("bi")) {
            CreateVCF();
            VcfProcedures.NumBiallelic(VCFreader, vcfwriter);

        }

        if (cmd.hasOption("MAF")) {
            frecrar = Double.parseDouble(cmd.getOptionValues("MAF")[0]);
            CreateVCF();
            VcfProcedures.Variantrare(VCFreader, frecrar, vcfwriter, variante);

        }

        if (cmd.hasOption("maxHet")) {
            minHet = Double.parseDouble(cmd.getOptionValues("maxHet")[0]);
            CreateVCF();
            VcfProcedures.MinHet(VCFreader, minGQ, vcfwriter, variante);
        }

        if (cmd.hasOption("nr")) {
            Size = Integer.parseInt(cmd.getOptionValues("nr")[0]);
            CreateVCF();
            VcfProcedures.BestQUALinKb(VCFreader, Size, vcfwriter);

        }

        //if (cmd.hasOption("split")){
        //    Split();
        //}
    }

    public static org.apache.commons.cli.Options CreateOptions() {

        org.apache.commons.cli.Options options = new org.apache.commons.cli.Options();

        Option input = OptionBuilder.withArgName("path to input VCF file")
                .hasArgs(1)
                .withDescription("Required, must be uncompressed. Will index it.")
                .create("input");

        Option output = OptionBuilder.withArgName("path to output VCF file")
                .hasArgs(1)
                .withDescription("Optional, by default prints to STDOUT.")
                .create("output");

        Option DP = OptionBuilder.withArgName("integer")
                .hasArgs(1)
                .withDescription("Minimum overall DP of each variant (row) in the input.")
                .create("DP");

        Option sDP = OptionBuilder.withArgName("integer")
                .hasArgs(1)
                .withDescription("Minimum DP of each sample.")
                .create("sDP");

        Option missing = OptionBuilder.withArgName("double")
                .hasArgs(1)
                .withDescription("Allowed % of missing data in each variant. Example: -missing 5")
                .create("missing");

        Option interval = OptionBuilder.withArgName("integer")
                .hasArgs(2)
                .withDescription("Select interval of coordinates, requires -chr . Example: -interval 1000 2000")
                .create("interval");

        Option chr = OptionBuilder.withArgName("String")
                .hasArgs(1)
                .withDescription("Select a chromosome of interest. Example: -chr Bd1")
                .create("chr");

        Option pos = OptionBuilder.withArgName("integer")
                .hasArgs(1)
                .withDescription("Select a precise position. Example: -pos 1234")
                .create("pos");

        Option sample = OptionBuilder.withArgName("String: sample names")
                .hasArgs()
                .withDescription("Select a subset of samples. Example: -sample sample1,sample2")
                .create("sample");

        Option call = OptionBuilder.withArgName("String: sample name")
                .hasArgs(1)
                .withDescription("Requires -pos & -sample. Example: -call sample1 -chr Bd1 -pos 1000")
                .create("call");

        Option bi = new Option("bi", "Select only biallelic genotypes.");

        Option MAF = OptionBuilder.withArgName("double")
                .hasArgs(1)
                .withDescription("Minimum allele frequency [0-1]. Example: -MAF 0.05")
                .create("MAF");

        Option maxHet = OptionBuilder.withArgName("double")
                .hasArgs(1)
                .withDescription("Maximum % of heterozigous samples. Example: -maxHet 90")
                .create("maxHet");

        Option nr = OptionBuilder.withArgName("integer in Kb")
                .hasArgs(1)
                .withDescription("Select top quality variant in selected window size, in Kb. Example: -nr 100")
                .create("nr");

        Option menu = new Option("menu", "Select step-by-step menu.");

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
        options.addOption(menu);
        //options.addOption(split);
        return options;
    }

    public static void main(String[] args) throws IOException, ParseException {

        org.apache.commons.cli.Options options = CreateOptions();

        CommandLineParser parser = new DefaultParser();

        cmd = parser.parse(options, args);

        if (args.length < 1) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", options);

            System.out.println("\nThis utility builds on HTSJDK and handles VCF versions supported there, currently v4.2.");
            System.out.println("Eduardo Candeal, with help from Carlos P Cantalapiedra and Bruno Contreras-Moreira\nEEAD-CSIC 2016");
        } else {
            VcfFilter vcffilter = new VcfFilter();
            vcffilter.ReadOptions();
        }
    }
}
