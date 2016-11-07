
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.EnumSet;
import java.util.HashSet;
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
    
    String pathin;
    public VCFFileReader VCFreader;
    
    public VcfController(){
        
    }
    
    
    
    public void ReadOptions(CommandLine cmd) throws IOException {

        String variante = "";

        // TODO: input option is required
//        if (cmd.hasOption("input")) {
        pathin = cmd.getOptionValues("input")[0];
        // TODO: check whether .idx file exists
        
//        }

//        if (pathin.isEmpty()) {
//            /* VCFreader = new VCFFileReader(System.in); // does not work as 
//            VCFFileReader constructor expects a java.io.File file, not a stream */
//        } else {
//            VCFreader = new VCFFileReader(new File(pathin));
//        }
        

        String pathout = ""; // System.out or specified path
        if (cmd.hasOption("output")) {
            pathout = cmd.getOptionValues("output")[0];
        }
        
        String idxFile = pathin + ".idx";
        File inputFile = new File(pathin);
        VcfUtils.createidx(idxFile, inputFile);
        
        VCFreader = new VCFFileReader(inputFile);
        VCFHeader header = VCFreader.getFileHeader();
        VariantContextWriter vcfwriter = VcfUtils.createVCF(header, pathout);

        if (cmd.hasOption("DP")) {
            int DPG = Integer.parseInt(cmd.getOptionValues("DP")[0]);
            VcfProcedures.MinDPGen(VCFreader, DPG, vcfwriter, variante);
        }
        if (cmd.hasOption("sDP")) {
            int DPS = Integer.parseInt(cmd.getOptionValues("sDP")[0]);
            VcfProcedures.MinDPSample(VCFreader, DPS, vcfwriter, variante);
        }

        if (cmd.hasOption("missing")) {
            double NData = Integer.parseInt(cmd.getOptionValues("missing")[0]);
            VcfProcedures.MissingData(VCFreader, NData, vcfwriter, variante);
        }

        String crom = "";
        if (cmd.hasOption("chr")) {
            crom = cmd.getOptionValues("chr")[0];
        }

        if (cmd.hasOption("interval")) {
            // TODO: must have "chr" option also
            crom = cmd.getOptionValues("chr")[0];
            int posfirstint = Integer.parseInt(cmd.getOptionValues("interval")[0]);
            int possecondint = Integer.parseInt(cmd.getOptionValues("interval")[1]);
            VcfProcedures.VariantInter(VCFreader, posfirstint, possecondint,
                    crom, vcfwriter, variante);
        }

        int position = -1;
        if (cmd.hasOption("pos")) {
            position = Integer.parseInt(cmd.getOptionValues("pos")[0]);
        }

        if (cmd.hasOption("sample")) {

            Set set = new HashSet();
            String[] samples2 = cmd.getOptionValues("sample")[0].split(",");
            System.out.println(samples2.length);
            for (int i = 0; i < samples2.length; i++) {
                System.out.println(samples2[i]);
                set.add(samples2[i]);
            }
            VcfProcedures.SelectGenotype(VCFreader, set, pathout);

        }

        if (cmd.hasOption("call")) {

            String sample = cmd.getOptionValues("call")[0];
            VcfProcedures.FindSamVar(VCFreader, position, crom,
                    sample, variante);
            System.out.println(variante);

        }
        if (cmd.hasOption("bi")) {
            VcfProcedures.NumBiallelic(VCFreader, vcfwriter);

        }

        if (cmd.hasOption("MAF")) {
            double frecrar = Double.parseDouble(cmd.getOptionValues("MAF")[0]);
            VcfProcedures.Variantrare(VCFreader, frecrar, vcfwriter, variante);

        }

        if (cmd.hasOption("maxHet")) {
            double minHet = Double.parseDouble(cmd.getOptionValues("maxHet")[0]);
            int minGQ = 0; // TODO: initialize
            VcfProcedures.MinHet(VCFreader, minGQ, vcfwriter, variante);
        }

        if (cmd.hasOption("nr")) {
            int Size = Integer.parseInt(cmd.getOptionValues("nr")[0]);
            VcfProcedures.BestQUALinKb(VCFreader, Size, vcfwriter);

        }

        //if (cmd.hasOption("split")){
        //    Split();
        //}
        
        
    }
}
