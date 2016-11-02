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
    public EnumSet<Options> DEFAULT_OPTIONS=
            EnumSet.of(Options.INDEX_ON_THE_FLY, 
            Options.ALLOW_MISSING_FIELDS_IN_HEADER, 
            Options.WRITE_FULL_FORMAT_FIELD);
    
    
    
    
    public void MinDPGen(){
               
        Iterator<VariantContext> iter=VCFreader.iterator();
        while(iter.hasNext()){
        VariantContext variant= iter.next();
        int DPV=variant.getAttributeAsInt("DP", 0);
        if(DPV>DPG){
            variante=variante+variant.toString()+"\n";
            vcfwriter.add(variant);
        }
        }
       
    }
    public void MinDPSample(){
 
       total=VCFreader.getFileHeader().getNGenotypeSamples();
        Iterator<VariantContext> iter=VCFreader.iterator();
        while(iter.hasNext()){
            VariantContext variant= iter.next();
            boolean minimo=false;
            for(int i=0; i<total;i++){
            int DP= variant.getGenotype(i).getDP();
               if(DPS>DP){
                   minimo=true;
               } 
            }
            if(minimo==false){
              variante=variante+variant.toString()+"\n";
              vcfwriter.add(variant);
            }
        }
     
     }   
     


     public void MinHet(){
         
         total=VCFreader.getFileHeader().getNGenotypeSamples();
         double maxhet=(double)minHet*(double)total/100;
          
         VCFreader.iterator().forEachRemaining(variantcontext -> {
             int numhet=0;
             
             for (int i=0; i<total; i++){
                 if(variantcontext.getGenotype(i).isHet()){
                 numhet=numhet+1;
                 }                
             }
             if (maxhet>numhet){
                 vcfwriter.add(variantcontext);
             }      
             
         });

         
        }
     
    public void Variantrare(){
           
           total=VCFreader.getFileHeader().getNGenotypeSamples();
            HashMap h=new HashMap();
           VCFreader.iterator().forEachRemaining(variantcontext -> {
            double totalalleles=0;
            h.clear();
               for (int i=0; i<total; i++){
                   for (Allele o: variantcontext.getGenotype(i).getAlleles()){
                   totalalleles=totalalleles+1;    
                   String b=o.getBaseString();
                   if(h.containsKey(b)){
                   h.put(b, (int)h.get(b)+1 );
                   }
                   else {
                   h.put(b, 1);
                   }
                   
                   
                   }

                 }
               boolean pass=true;
               Iterator iter=h.keySet().iterator();
               while(iter.hasNext()){        
               String c=(String)iter.next();
               double frec=(int)h.get(c)/totalalleles;
               if(frec<frecrar){
                   pass=false;
               }
               }
               if(pass==true){
               vcfwriter.add(variantcontext);
               }
             

             }); 
         
        } 

    public void MissingData(){
        
       
       total=VCFreader.getFileHeader().getNGenotypeSamples();
       NData=total*NData/100;


       VCFreader.iterator().forEachRemaining(variantcontext -> {
       int num=0;
       for(int i=0; i<total;i++){
      
        int DP2= variantcontext.getGenotype(i).getDP();
           if(DP2==0){
               num=num+1;
           } 
           
        }           
        if(num<NData){
            vcfwriter.add(variantcontext);
        }
        });
       
       
       }
    
    public void NumBiallelic(){
           
           VCFreader.iterator().forEachRemaining(variantcontext -> {
    
            if(variantcontext.isSNP() && variantcontext.isBiallelic()){
                vcfwriter.add(variantcontext);
            
           } 
        });
       }
    
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
    
    public void FindSamVar(){
       
       Iterator<VariantContext> iter=VCFreader.iterator();
       boolean find=false;

        while(iter.hasNext()){
            VariantContext variant= iter.next();
            int num=variant.getEnd();
            String name=variant.getContig();
            if(num==position && name.equals(crom)){
                variante=variant.getGenotype(sample).toString();
                find=true;
                break;
             
            }    
        }
        if(find==false){
        System.out.println("your sample doesn't exit");
        }

    }
          
    
    
    public void SelectGenotype() throws IOException {
       
         VCFHeader header= new VCFHeader(VCFreader.getFileHeader().getMetaDataInInputOrder(), set);
         FileOutputStream outputstream= new FileOutputStream(new File(pathout));
         vcfwriter=VariantContextWriterFactory.createVcf(new File(pathout), 
                outputstream, VCFreader.getFileHeader().getSequenceDictionary(),DEFAULT_OPTIONS);
     
          vcfwriter.writeHeader(header); 
            
        Iterator<VariantContext> iter= VCFreader.iterator();
        while(iter.hasNext()){
        VariantContext variant=iter.next();
         VariantContext vc2=variant.subContextFromSamples(set);
                vcfwriter.add(vc2);
        }
    
    }

    public void VariantInter(){

       
        Iterator<VariantContext> iter=VCFreader.iterator();
        
        System.err.println("VariantInter "+posfirstint+"-"+possecondint);
        
        while(iter.hasNext()){
            VariantContext variant= iter.next();
            String id= variant.getContig();
            
            if (variant.isSNP() &&
                variant.getEnd()>=posfirstint &&
                variant.getEnd()<=possecondint &&
                id.equals(crom)){
                
                vcfwriter.add(variant);
                variante=variante+variant.toString()+"\n";
                
            }
            
            // ASSERT: VCF INPUT FILE IS SORTED BY CHR AND POSITION
            if(variant.getEnd()>possecondint && id.equals(crom)){
                break;
            }
        }

    }
    
    public void Split(){
        Split splits= new Split();
        splits.ifile=pathin;
        splits.init();
    }
    
    public void BestQUALinKb(){

        Iterator<VariantContext> ita= VCFreader.iterator();  
        ita.next();
        VariantContext variantref= ita.next();

        int posref= variantref.getEnd()+Size;
        String cromo=variantref.getContig();
        double qual=variantref.getPhredScaledQual();
        int posvariant=0;
        int i=0;
        VariantContext max=variantref;
        while(ita.hasNext()){
            
            VariantContext variant =ita.next();
            posvariant=variant.getEnd();
            
                if(cromo.equals(variant.getContig())==false){
                    vcfwriter.add(max);
                    cromo=variant.getContig();
                    posref=variant.getEnd()+Size;
                    qual=variant.getPhredScaledQual();
                    max=variant;
                
                }
            
                if (posvariant<posref && variant.getContig().equals(cromo) && variant.getPhredScaledQual()>qual){
                    qual=variant.getPhredScaledQual();
                    max=variant;
                    
                }

                if(posvariant>posref && variant.getContig().equals(cromo)){
                    vcfwriter.add(max);
                    posref=variant.getEnd()+Size;
                    qual=variant.getPhredScaledQual();
                    max=variant;
                
                }

             i++;

        }
    
    }
    
    public void MinGQkb(){

        total=VCFreader.getFileHeader().getNGenotypeSamples();
        ArrayList<Integer> array= new ArrayList();

        Iterator<VariantContext> ita= VCFreader.iterator();  
        ita.next();
        VariantContext variantref= ita.next();

        int posref= variantref.getEnd()+Size*1000;
        String cromo=variantref.getContig();
        int num=0;
        int posvariant=0;
        while(ita.hasNext()){
            
            VariantContext variant =ita.next();
            boolean find=false;
            posvariant=variant.getEnd();
                if (posvariant<posref && variant.getContig()==cromo){
                    for (int i=0; i<total;i++){
                        if(variant.getGenotype(i).getGQ()<minGQ){
                            find=true;      
                        }
                        
                    } 
                    if(find==false){
                    num+=1;
                    }
                              
                }
                else {
                posref=variant.getEnd()+Size*1000;
                cromo=variant.getContig();
                array.add(num);
                num=0;
                
                }
        }
        variante=""+array.size();
  
    
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
        } finally {
            if (stream != null) {
                stream.close();
            }
        }
    }
    
    
    public void createidx() {
        int binSize = 16000;
        File inputFile = new File(pathin);
        VCFCodec codec = new VCFCodec();
        
        String idxFile=pathin+".idx";
      
        AbstractIndex idx= IndexFactory.createLinearIndex(inputFile, codec, binSize); 
       
        //try{
        writeTribbleIndex(idx,idxFile);      
        //} catch (Exception e){}
    }
    
    public void ReadOptions() throws IOException {
       
        if (cmd.hasOption("input")){
            pathin=cmd.getOptionValues("input")[0];
            createidx();
        }
        
        if (pathin.isEmpty()){  
            /* VCFreader = new VCFFileReader(System.in); // does not work as 
            VCFFileReader constructor expects a java.io.File file, not a stream */
        }
        else {
            VCFreader = new VCFFileReader(new File(pathin));
        }
           

        if (cmd.hasOption("output")){
              pathout=cmd.getOptionValues("output")[0];
        } 
        
        if(cmd.hasOption("DP")) {  
                DPG=Integer.parseInt(cmd.getOptionValues("DP")[0]);
                CreateVCF();
                MinDPGen();  
            }
        if (cmd.hasOption("sDP")) {     
                DPS=Integer.parseInt(cmd.getOptionValues("sDP")[0]);
                CreateVCF();
                MinDPSample(); 
        }
        
        if (cmd.hasOption("missing")) {     
                NData=Integer.parseInt(cmd.getOptionValues("missing")[0]);
                CreateVCF();
                MissingData();
        }
        
        if (cmd.hasOption("chr")) {
           crom=cmd.getOptionValues("chr")[0];
        }
        
        if (cmd.hasOption("interval")){
                posfirstint=Integer.parseInt(cmd.getOptionValues("interval")[0]);
                possecondint=Integer.parseInt(cmd.getOptionValues("interval")[1]);
                CreateVCF();
                VariantInter();
        }
        
        if(cmd.hasOption("pos")){
            position=Integer.parseInt(cmd.getOptionValues("pos")[0]);
        }
        
        if (cmd.hasOption("sample")){
                        
                set= new HashSet();
                String [] samples2=cmd.getOptionValues("sample")[0].split(",");
                System.out.println(samples2.length);
                for (int i=0; i<samples2.length; i++){
                    System.out.println(samples2[i]);
                set.add(samples2[i]);
                }
                SelectGenotype();
        
        }
        
        if (cmd.hasOption("call")){
        
                sample=cmd.getOptionValues("call")[0];
                FindSamVar();
                System.out.println(variante);
        
        }
        if (cmd.hasOption("bi")){
                CreateVCF();
                NumBiallelic(); 
        
        }
        
        if (cmd.hasOption("MAF")){
                frecrar=Double.parseDouble(cmd.getOptionValues("MAF")[0]);
                CreateVCF();
                Variantrare();
        
        }
        
        if (cmd.hasOption("maxHet")){
                minHet=Double.parseDouble(cmd.getOptionValues("maxHet")[0]);
                CreateVCF();
                MinHet();
        }
        
        if (cmd.hasOption("nr")){
               Size=Integer.parseInt(cmd.getOptionValues("nr")[0]);
                CreateVCF();
                BestQUALinKb();
        
        }
        
        if(cmd.hasOption("menu")){
            menu();
        
        }
        
        //if (cmd.hasOption("split")){
        //    Split();
        //}
        

        
    }
    
    public static org.apache.commons.cli.Options CreateOptions(){
        
        org.apache.commons.cli.Options options = new org.apache.commons.cli.Options();
        
        Option input = OptionBuilder.withArgName( "path to input VCF file" )                    
                                .hasArgs(1)
                                .withDescription("Required, must be uncompressed. Will index it." )
                                .create( "input" );
        
        Option output = OptionBuilder.withArgName( "path to output VCF file" )                    
                                .hasArgs(1)
                                .withDescription("Optional, by default prints to STDOUT.")
                                .create( "output" );        
        
        Option DP = OptionBuilder.withArgName( "integer" )
                                .hasArgs(1)
                                .withDescription("Minimum overall DP of each variant (row) in the input.")
                                .create( "DP" );
        
        Option sDP = OptionBuilder.withArgName( "integer" )
                                .hasArgs(1)
                                .withDescription("Minimum DP of each sample.")
                                .create( "sDP" );
        
        Option missing = OptionBuilder.withArgName( "double" )
                                .hasArgs(1)
                                .withDescription("Allowed % of missing data in each variant. Example: -missing 5")
                                .create( "missing" );
        
        Option interval = OptionBuilder.withArgName( "integer" )
                                .hasArgs(2)
                                .withDescription("Select interval of coordinates, requires -chr . Example: -interval 1000 2000")
                                .create( "interval" );
        
        Option chr = OptionBuilder.withArgName( "String" )
                                .hasArgs(1)
                                .withDescription( "Select a chromosome of interest. Example: -chr Bd1" )
                                .create( "chr" );
        
        Option pos = OptionBuilder.withArgName( "integer" )
                                .hasArgs(1)
                                .withDescription("Select a precise position. Example: -pos 1234")
                                .create( "pos" );
        
        Option sample = OptionBuilder.withArgName( "String: sample names" )
                                .hasArgs()
                                .withDescription("Select a subset of samples. Example: -sample sample1,sample2" )
                                .create( "sample" );
        
        Option call = OptionBuilder.withArgName( "String: sample name" )
                                .hasArgs(1)
                                .withDescription("Requires -pos & -sample. Example: -call sample1 -chr Bd1 -pos 1000" )
                                .create( "call" );
        
        Option bi = new Option( "bi", "Select only biallelic genotypes." );
        
        Option MAF = OptionBuilder.withArgName("double")
                                .hasArgs(1)
                                .withDescription("Minimum allele frequency [0-1]. Example: -MAF 0.05" )
                                .create( "MAF" );
        
        Option maxHet = OptionBuilder.withArgName("double")
                                .hasArgs(1)
                                .withDescription("Maximum % of heterozigous samples. Example: -maxHet 90" )
                                .create( "maxHet" );
        
        Option nr = OptionBuilder.withArgName("integer in Kb")
                                .hasArgs(1)
                                .withDescription("Select top quality variant in selected window size, in Kb. Example: -nr 100" )
                                .create( "nr" );
        
        Option menu = new Option( "menu", "Select step-by-step menu." );
        
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
        
        if(args.length<1) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", options);
            
            System.out.println("\nThis utility builds on HTSJDK and handles VCF versions supported there, currently v4.2.");
            System.out.println("Eduardo Candeal, with help from Carlos P Cantalapiedra and Bruno Contreras-Moreira\nEEAD-CSIC 2016");
        } else {
            VcfReader vcffile = new VcfReader();
            vcffile.ReadOptions();
        }
    }
}
