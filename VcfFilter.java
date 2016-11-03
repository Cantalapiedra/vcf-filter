
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
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
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
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();

        if (pathout.isEmpty()) {
            builder.setOutputFile((File)null);
            builder.setOutputStream(System.out);
//            vcfwriter = VariantContextWriterFactory.createVcf(null, System.out,
//                    VCFreader.getFileHeader().getSequenceDictionary(), DEFAULT_OPTIONS);
        } else {
            FileOutputStream outputstream = new FileOutputStream(new File(pathout));

            builder.setOutputFile(new File(pathout));
            builder.setOutputStream(outputstream);
            
//            vcfwriter = VariantContextWriterFactory.createVcf(new File(pathout), outputstream,
//                    VCFreader.getFileHeader().getSequenceDictionary(), DEFAULT_OPTIONS);
        }
        
        builder.setReferenceDictionary(VCFreader.getFileHeader().getSequenceDictionary());
        builder.setOptions(DEFAULT_OPTIONS);
        vcfwriter = builder.build();
        
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
