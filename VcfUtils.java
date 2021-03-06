
import htsjdk.tribble.index.AbstractIndex;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.OutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.EnumSet;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author CPCantalapiedra 2016
 */
public class VcfUtils {
    
    public static EnumSet<Options> DEFAULT_OPTIONS
            = EnumSet.of(Options.INDEX_ON_THE_FLY,
                    Options.ALLOW_MISSING_FIELDS_IN_HEADER,
                    Options.WRITE_FULL_FORMAT_FIELD);
    
    public static EnumSet<Options> STREAM_OPTIONS
            = EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER,
                    Options.WRITE_FULL_FORMAT_FIELD);
    
    public static VariantContextWriter createVCF(VCFHeader header, String pathout) throws IOException {
        
        File outfile = null;
        OutputStream outputstream = System.out;
        
        if (!pathout.isEmpty()) {
            outfile = new File(pathout);
            outputstream = new FileOutputStream(outfile);
        }
        
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(outfile);
        builder.setOutputStream(outputstream);

        builder.setReferenceDictionary(header.getSequenceDictionary());
        builder.setOptions(STREAM_OPTIONS);
        VariantContextWriter vcfwriter = builder.build();

        return vcfwriter;
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

    public static void createidx(String idxFilename, File vcfFile) throws IOException {
        int binSize = 16000;
        VCFCodec codec = new VCFCodec();

        AbstractIndex idx = IndexFactory.createLinearIndex(vcfFile, codec, binSize);

        //try{
        writeTribbleIndex(idx, idxFilename);
        //} catch (Exception e){}
    }
}
