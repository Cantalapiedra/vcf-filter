// CPCantalapiedra 2016

import java.util.Iterator;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

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

import java.util.Scanner;
import java.util.Set;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;

public class VcfProcedures {
    
    public void MinDPGen(VCFFileReader VCFreader, int DPG, VariantContextWriter vcfwriter){
        retValue = false;
        
        Iterator<VariantContext> iter = VCFreader.iterator();
        VariantContext variant;
        int DPV;
        while(iter.hasNext()){
            variant = iter.next();
            DPV = variant.getAttributeAsInt("DP", 0);
            
            if(DPV > DPG){
                retValue = true;
                variante = variante+variant.toString()+"\n";
                vcfwriter.add(variant);
            }
        }
        
        return retValue;
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
    
}