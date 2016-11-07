
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import org.apache.commons.cli.CommandLine;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author CPCantalapiedra 2016
 */
public class VcfController {
    
    CommandLine cmd;
    
    String pathin;
    String pathout = ""; // STDOUT if empty
    
    VCFFileReader vcfReader;
    VCFHeader header;
    VariantContextWriter vcfwriter;
    
    // Options
    String chr = "";
    int intervalStart = -1;
    int intervalEnd = -1;
    int pos = -1;
    Set samplesSet;
    int variantDP = -1;
    int DPS;
    double missPercen = -1.0;
    String sampleCall;
    boolean biallelic = false;
    double frecrar;
    double minHet;
    int minGQ;
    int Size;
    
    public VcfController(CommandLine cmd){
        this.cmd = cmd;
    }
    
    public void run() throws IOException, Exception {
        readOptions();
        createFiles(pathin, pathout);
        ArrayList<VariantContext> variants = process();
        output(variants);
    }
    
    public void readOptions() throws Exception {
        
        /* Input and output files
        */
        pathin = cmd.getOptionValue("input");
        System.err.println("Input file: "+pathin);
        
        if (cmd.hasOption("output")) {
            pathout = cmd.getOptionValue("output");
            System.err.println("Output file: "+pathout);
        } else {
            System.err.println("Output file: STDOUT");
        }
        
        /***** Positional filtering 
        */
        // Chromosome
        if (cmd.hasOption("chr")) {
            chr = cmd.getOptionValue("chr");
            System.err.println("Chromosome: "+chr);
        }
        
        // Position
        if (cmd.hasOption("pos")) {
            if (!cmd.hasOption("chr")) // chr is required
                throw new Exception("chr option is required for position.");
            
            pos = Integer.parseInt(cmd.getOptionValues("pos")[0]);
            System.err.println("Position: "+chr);
        }
        
        // Range of positions
        if (cmd.hasOption("r")) {
            if (!cmd.hasOption("chr")) // chr is required
                throw new Exception("chr option is required for range.");
            
            String[] interval = cmd.getOptionValues("r");
            intervalStart = Integer.parseInt(interval[0]);
            intervalEnd = Integer.parseInt(interval[1]);
            System.err.println("Range: "+intervalStart+"-"+intervalEnd);
        }
        
        /****** By variant data
         */
        // Depth of variant
        if (cmd.hasOption("dp"))
            variantDP = Integer.parseInt(cmd.getOptionValues("dp")[0]);
        
        // Missing data
        if (cmd.hasOption("m"))
            missPercen = Integer.parseInt(cmd.getOptionValues("missing")[0]);
        
        
        if (cmd.hasOption("call")) {
            if (!cmd.hasOption("pos")) // pos is required
            sampleCall = cmd.getOptionValues("call")[0];
        }
        
        // Specific samples
        if (cmd.hasOption("sample")) {

            samplesSet = new HashSet();
            String[] samples2 = cmd.getOptionValues("sample")[0].split(",");
            System.err.println("Number of samples specified: "+samples2.length);
            for (int i = 0; i < samples2.length; i++) {
                System.err.println(samples2[i]);
                samplesSet.add(samples2[i]);
            }
        }
        
        if (cmd.hasOption("sDP"))
            DPS = Integer.parseInt(cmd.getOptionValues("sDP")[0]);

        
        
        if (cmd.hasOption("bi"))
            biallelic = true;

        if (cmd.hasOption("MAF"))
            frecrar = Double.parseDouble(cmd.getOptionValues("MAF")[0]);

        if (cmd.hasOption("maxHet")) {
            minHet = Double.parseDouble(cmd.getOptionValues("maxHet")[0]);
            minGQ = 0; // TODO: initialize
        }

        if (cmd.hasOption("nr")) {
            Size = Integer.parseInt(cmd.getOptionValues("nr")[0]);
        }

        //if (cmd.hasOption("split")){
        //    Split();
        //}
    }
    
    public void createFiles(String pathin, String pathout) throws IOException {
        // TODO: check whether .idx file exists
        String idxFile = pathin + ".idx";
        File inputFile = new File(pathin);
        VcfUtils.createidx(idxFile, inputFile);
        
        vcfReader = new VCFFileReader(inputFile);
        header = vcfReader.getFileHeader();
        vcfwriter = VcfUtils.createVCF(header, pathout);
    }
    
    public ArrayList<VariantContext> process() throws IOException {
        
        String variante = "";
        
        // List to store the variants to be output
        ArrayList<VariantContext> variants = new ArrayList();
        
        // Read variants of input VCF/BCF file
        Iterator<VariantContext> iter = vcfReader.iterator();

        VariantContext variant;
        boolean chrFound = false;
        while (iter.hasNext()) {
            variant = iter.next();
            
            // Chromosome
            if (!chr.isEmpty()){
                if (!VcfProcedures.variantChr(variant, chr))
                    continue;
            
                // ASSERT: VCF INPUT FILE IS SORTED BY CHR
                if (variant.getContig().equals(chr))
                    chrFound = true;
                if ((!variant.getContig().equals(chr)) && chrFound == true)
                        break;
            }
            
            // Specific position
            if (pos != -1){
                if (!VcfProcedures.variantPos(variant, pos))
                    continue;
            }
            
            // Interval or range of positions
            if (intervalStart != -1 && intervalEnd != -1){
                if (!VcfProcedures.variantRange(variant, 
                        intervalStart, intervalEnd)) {
                    // ASSERT: VCF INPUT FILE IS SORTED BY CHR AND POSITION
                    if (variant.getStart() > intervalEnd && 
                            variant.getContig().equals(chr))
                        break;
                    else
                        continue;
                }
            }
            
            // Filter by general DP of variant
            if (variantDP != -1){
                if (!VcfProcedures.minVariantDP(variant, variantDP))
                    continue;
            }
            
            // Missing data
            if (missPercen != -1) {
                if (!VcfProcedures.missingData(variant, missPercen, header.getNGenotypeSamples()))
                    continue;
            }
            
            variants.add(variant);
        }
        
        
        
        if (cmd.hasOption("sDP")) {
            VcfProcedures.MinDPSample(vcfReader, DPS, vcfwriter, variante);
        }

        
        
        

        if (cmd.hasOption("sample")) {

            Set set = new HashSet();
            String[] samples2 = cmd.getOptionValues("sample")[0].split(",");
            System.out.println(samples2.length);
            for (int i = 0; i < samples2.length; i++) {
                System.out.println(samples2[i]);
                set.add(samples2[i]);
            }
            VcfProcedures.SelectGenotype(vcfReader, set, pathout);

        }

        if (cmd.hasOption("call")) {
            VcfProcedures.FindSamVar(vcfReader, pos, chr,
                    sampleCall, variante);
            System.out.println(variante);

        }
        if (cmd.hasOption("bi")) {
            VcfProcedures.NumBiallelic(vcfReader, vcfwriter);

        }

        if (cmd.hasOption("MAF")) {
            VcfProcedures.Variantrare(vcfReader, frecrar, vcfwriter, variante);

        }

        if (cmd.hasOption("maxHet")) {
            VcfProcedures.MinHet(vcfReader, minGQ, vcfwriter, variante);
        }

        if (cmd.hasOption("nr")) {
            VcfProcedures.BestQUALinKb(vcfReader, Size, vcfwriter);
        }
        
        return variants;
    }
    
    public void output(ArrayList<VariantContext> variants) {
        vcfwriter.writeHeader(header);
        variants.forEach(variant -> vcfwriter.add(variant));
    }
}
