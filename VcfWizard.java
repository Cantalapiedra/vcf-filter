/**
 *
 * @author Eduardo Candeal 2016
 */

import htsjdk.tribble.index.AbstractIndex;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.io.IOException;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Set;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;

public class VcfReader {

    public int num;
    public String ifile;
    public String pathin="";
    public String pathout="";

    public VCFFileReader VCFreader;
    
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
    public EnumSet<Options> DEFAULT_OPTIONS=
            EnumSet.of(Options.INDEX_ON_THE_FLY, 
            Options.ALLOW_MISSING_FIELDS_IN_HEADER, 
            Options.WRITE_FULL_FORMAT_FIELD);
    
    public void CreateVCF()throws IOException {
        VCFreader =new VCFFileReader(new File(pathin));
        if(pathout.isEmpty()){
                vcfwriter=VariantContextWriterFactory.createVcf(null, System.out,
                VCFreader.getFileHeader().getSequenceDictionary(),DEFAULT_OPTIONS);
        }
        else {
                FileOutputStream outputstream= new FileOutputStream(new File(pathout));
                vcfwriter=VariantContextWriterFactory.createVcf(new File(pathout), outputstream,
                VCFreader.getFileHeader().getSequenceDictionary(),DEFAULT_OPTIONS);

        }

        vcfwriter.writeHeader(VCFreader.getFileHeader());
    }
    
    
    
    public void menu() throws IOException {
        
        Scanner sc= new Scanner(System.in);   
        System.out.println("First, type the path to your INPUT file");             
        pathin=sc.next();
        System.out.println("Now, type the path to your OUTPUT file");
        pathout=sc.next();
        int opt;
        do {
        System.out.println("-----------------------------");
        System.out.println("Options Menu in VCF File");
        System.out.println("-----------------------------");
        System.out.println("Option 1- Filter by minimun General DP in VariantContext");
        System.out.println("Option 2- Filter by minimun  DP in each Genotype");
        System.out.println("Option 3- Filter by maximun of missing data of genotypes");
        System.out.println("Option 4- Get a sample in a specific position");
        System.out.println("Option 5- Obtain number of SNPS with minimun  sample GQ in each Genotype in a specific distance in kb");
        System.out.println("Option 6- Select SNPs in an interval in a concrete cromosome");
        System.out.println("Option 7- Filter by Biallelic SNPs");
        System.out.println("Option 8- Filter by minimum allele frequency");
        System.out.println("Option 9- Filter by % maximun of Heterozigotes samples");
        //System.out.println("Option 10- Split VCF File in to SNPs and Indels");
        System.out.println("Option 11- Select samples and generate new vcf file");
        System.out.println("Option 12- Select the best quality variant in your desired distance");
        System.out.println("Press 0 to leave");
        System.out.println("----------------------------------------------------");   
        System.out.println("Choose the option number");
        opt=sc.nextInt();
        
        switch (opt) {

            case 1:
                System.out.println("Type the minimum General DP");
                DPG=sc.nextInt();
                CreateVCF();
                MinDPGen();

            break;
            
            case 2:
                System.out.println("Type the minimum sample DP");
                DPS=sc.nextInt();
                CreateVCF();
                MinDPSample();
            
            break;
                
            case 3:
                System.out.println("Type the % of maximum missing data per variantcontext");
                NData=sc.nextInt();
                CreateVCF();
                MissingData();
            
            break;
            
            case 4:
                System.out.println("Type the sample to find");
                sample=sc.next();
                System.out.println("Type the specific position");
                position=sc.nextInt();
                System.out.println("Type the specific chromosome");
                crom=sc.next();
                FindSamVar();
                System.out.println(variante);
            
            break;
                        
            case 5:
                System.out.println("Type the minimun sample GQ");
                minGQ=sc.nextInt();
                System.out.println("Type the distance in Kb");
                Size=sc.nextInt();
                CreateVCF();
                MinGQkb();
                System.out.println(variante);
            
            break;
            
            case 6:
                System.out.println("Introduce the first position");
                posfirstint=sc.nextInt();
                System.out.println("Introduce the last position");
                possecondint=sc.nextInt();
                System.out.println("Introduce the specific cromosome");
                crom=sc.next();
                CreateVCF();
                VariantInter();
            
            break;
            
            case 7:
                CreateVCF();
                NumBiallelic();       
            
            break;
            
            case 8:
                System.out.println("Introduce the number of minimum allele frequency. Use , instead of .");
                frecrar=sc.nextDouble();
                CreateVCF();
                Variantrare();
            
            break;
            
            case 9:
                System.out.println("Introduce % of maximum heterozigotes samples");
                minHet=sc.nextDouble();
                CreateVCF();
                MinHet();
            
            break;
            
            //case 10:
            //    Split();
            //break;
            
            case 11:
                System.out.println("Write selected samples separated by ;");
                String samples1=sc.next();
                set= new HashSet();
                String [] samples2=samples1.split(";");
                System.out.println(samples2.length);
                for (int i=0; i<samples2.length; i++){
                    System.out.println(samples2[i]);
                set.add(samples2[i]);
                }
                SelectGenotype();
            break;
            
            case 12:
                System.out.println("Write your desired distance");
                Size=sc.nextInt();
                CreateVCF();
                BestQUALinKb();
            break;
            
            case 0:
            break;    
            
            default:
                System.out.println("Introduce correct option");    
            break;
   
        }         
        } while (opt!=0);
    }
    
    public static void main(String[] args) throws IOException, ParseException {
        
        if(args.length<1) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", options);
            
            System.out.println("\nThis utility builds on HTSJDK and handles VCF versions supported there, currently v4.2.");
            System.out.println("Eduardo Candeal, with help from Carlos P Cantalapiedra and Bruno Contreras-Moreira\nEEAD-CSIC 2016");
        } else {
            VcfReader vcffile = new VcfReader();
            vcffile.menu();
        }
    }
}
