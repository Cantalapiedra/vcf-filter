
import java.util.Iterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.io.IOException;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Set;

/**
 *
 * @author CPCantalapiedra 2016
 */
public class VcfProcedures {

    public static boolean MinDPGen(VCFFileReader VCFreader, int DPG,
            VariantContextWriter vcfwriter, String variante) {

        boolean retValue = false;

        Iterator<VariantContext> iter = VCFreader.iterator();
        VariantContext variant;
        int DPV;
        while (iter.hasNext()) {
            variant = iter.next();
            DPV = variant.getAttributeAsInt("DP", 0);

            if (DPV > DPG) {
                retValue = true;
                variante = variante + variant.toString() + "\n";
                vcfwriter.add(variant);
            }
        }

        return retValue;
    }

    public static void MinDPSample(VCFFileReader VCFreader, int DPS,
            VariantContextWriter vcfwriter, String variante) {

        int total = VCFreader.getFileHeader().getNGenotypeSamples();
        Iterator<VariantContext> iter = VCFreader.iterator();

        VariantContext variant;
        boolean minimo;
        int DP;

        while (iter.hasNext()) {
            variant = iter.next();
            minimo = false;

            for (int i = 0; i < total; i++) {
                DP = variant.getGenotype(i).getDP();
                if (DPS > DP) {
                    minimo = true;
                }
            }
            if (minimo == false) {
                variante = variante + variant.toString() + "\n";
                vcfwriter.add(variant);
            }
        }
    }

    public static void MinHet(VCFFileReader VCFreader, int minHet,
            VariantContextWriter vcfwriter, String variante) {

        int total = VCFreader.getFileHeader().getNGenotypeSamples();
        double maxhet = (double) minHet * (double) total / 100;

        VCFreader.iterator().forEachRemaining(variantcontext -> {
            int numhet = 0;

            for (int i = 0; i < total; i++) {
                if (variantcontext.getGenotype(i).isHet()) {
                    numhet = numhet + 1;
                }
            }
            if (maxhet > numhet) {
                vcfwriter.add(variantcontext);
            }

        }
        );
    }

    public static void Variantrare(VCFFileReader VCFreader, double frecrar,
            VariantContextWriter vcfwriter, String variante) {

        int total = VCFreader.getFileHeader().getNGenotypeSamples();
        HashMap h = new HashMap();

        VCFreader.iterator().forEachRemaining(variantcontext
                -> {
            double totalalleles = 0;
            h.clear();
            for (int i = 0; i < total; i++) {
                for (Allele o : variantcontext.getGenotype(i).getAlleles()) {
                    totalalleles = totalalleles + 1;
                    String b = o.getBaseString();
                    if (h.containsKey(b)) {
                        h.put(b, (int) h.get(b) + 1);
                    } else {
                        h.put(b, 1);
                    }

                }

            }
            boolean pass = true;
            Iterator iter = h.keySet().iterator();
            while (iter.hasNext()) {
                String c = (String) iter.next();
                double frec = (int) h.get(c) / totalalleles;
                if (frec < frecrar) {
                    pass = false;
                }
            }
            if (pass == true) {
                vcfwriter.add(variantcontext);
            }

        }
        );

    }

    public static void MissingData(VCFFileReader VCFreader, double NData,
            VariantContextWriter vcfwriter, String variante) {

        int total = VCFreader.getFileHeader().getNGenotypeSamples();
        double NData_per = total * NData / 100;

        VCFreader.iterator().forEachRemaining(variantcontext -> {
            int num = 0;
            for (int i = 0; i < total; i++) {

                int DP2 = variantcontext.getGenotype(i).getDP();
                if (DP2 == 0) {
                    num = num + 1;
                }
            }
            if (num < NData_per) {
                vcfwriter.add(variantcontext);
            }
        });

    }

    public static void NumBiallelic(VCFFileReader VCFreader, VariantContextWriter vcfwriter) {

        VCFreader.iterator().forEachRemaining(variantcontext -> {

            if (variantcontext.isSNP() && variantcontext.isBiallelic()) {
                vcfwriter.add(variantcontext);

            }
        }
        );
    }

    public static void FindSamVar(VCFFileReader VCFreader, int position, 
            String crom, String sample, String variante) {

        Iterator<VariantContext> iter = VCFreader.iterator();
        boolean find = false;

        while (iter.hasNext()) {
            VariantContext variant = iter.next();
            int num = variant.getEnd();
            String name = variant.getContig();
            if (num == position && name.equals(crom)) {
                variante = variant.getGenotype(sample).toString();
                find = true;
                break;

            }
        }
        if (find == false) {
            System.out.println("your sample doesn't exit");
        }

    }

    public static void SelectGenotype(VCFFileReader VCFreader, Set set, String pathout) throws IOException {

        VCFHeader header = new VCFHeader(VCFreader.getFileHeader().getMetaDataInInputOrder(), set);
        
        VariantContextWriter vcfwriter = VcfUtils.createVCF(header, pathout);
        
////        FileOutputStream outputstream = new FileOutputStream(new File(pathout));
////        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
////        
////        builder.setOutputFile(new File(pathout));
////        builder.setOutputStream(outputstream);
////        
////        builder.setReferenceDictionary(header.getSequenceDictionary());
////        builder.setOptions(DEFAULT_OPTIONS);
////        VariantContextWriter vcfwriter = builder.build();
//        
////        vcfwriter = VariantContextWriterFactory.createVcf(new File(pathout),
////                outputstream, VCFreader.getFileHeader().getSequenceDictionary(), DEFAULT_OPTIONS);
//
//        vcfwriter.writeHeader(header);

        Iterator<VariantContext> iter = VCFreader.iterator();
        while (iter.hasNext()) {
            VariantContext variant = iter.next();
            VariantContext vc2 = variant.subContextFromSamples(set);
            vcfwriter.add(vc2);
        }
    }

    public static void VariantInter(VCFFileReader VCFreader, int posfirstint,
            int possecondint, String crom, VariantContextWriter vcfwriter,
            String variante) {

        Iterator<VariantContext> iter = VCFreader.iterator();

        System.err.println("VariantInter " + posfirstint + "-" + possecondint);

        while (iter.hasNext()) {
            VariantContext variant = iter.next();
            String id = variant.getContig();

            if (variant.isSNP()
                    && variant.getEnd() >= posfirstint
                    && variant.getEnd() <= possecondint
                    && id.equals(crom)) {

                vcfwriter.add(variant);
                variante = variante + variant.toString() + "\n";

            }

            // ASSERT: VCF INPUT FILE IS SORTED BY CHR AND POSITION
            if (variant.getEnd() > possecondint && id.equals(crom)) {
                break;
            }
        }

    }

    public static void Split(String pathin) {
        Split splits = new Split();
        splits.ifile = pathin;
        splits.init();
    }

    public static void BestQUALinKb(VCFFileReader VCFreader, int Size,
            VariantContextWriter vcfwriter) {

        Iterator<VariantContext> ita = VCFreader.iterator();
        ita.next();
        VariantContext variantref = ita.next();

        int posref = variantref.getEnd() + Size;
        String cromo = variantref.getContig();
        double qual = variantref.getPhredScaledQual();
        int posvariant = 0;
        int i = 0;
        VariantContext max = variantref;
        while (ita.hasNext()) {

            VariantContext variant = ita.next();
            posvariant = variant.getEnd();

            if (cromo.equals(variant.getContig()) == false) {
                vcfwriter.add(max);
                cromo = variant.getContig();
                posref = variant.getEnd() + Size;
                qual = variant.getPhredScaledQual();
                max = variant;

            }

            if (posvariant < posref && variant.getContig().equals(cromo) && variant.getPhredScaledQual() > qual) {
                qual = variant.getPhredScaledQual();
                max = variant;

            }

            if (posvariant > posref && variant.getContig().equals(cromo)) {
                vcfwriter.add(max);
                posref = variant.getEnd() + Size;
                qual = variant.getPhredScaledQual();
                max = variant;

            }

            i++;

        }

    }

    public static void MinGQkb(VCFFileReader VCFreader, int Size, int minGQ,
            String variante) {

        int total = VCFreader.getFileHeader().getNGenotypeSamples();
        ArrayList<Integer> array = new ArrayList();

        Iterator<VariantContext> ita = VCFreader.iterator();
        ita.next();
        VariantContext variantref = ita.next();

        int posref = variantref.getEnd() + Size * 1000;
        String cromo = variantref.getContig();
        int num = 0;
        int posvariant = 0;
        while (ita.hasNext()) {

            VariantContext variant = ita.next();
            boolean find = false;
            posvariant = variant.getEnd();
            if (posvariant < posref && variant.getContig() == cromo) {
                for (int i = 0; i < total; i++) {
                    if (variant.getGenotype(i).getGQ() < minGQ) {
                        find = true;
                    }

                }
                if (find == false) {
                    num += 1;
                }

            } else {
                posref = variant.getEnd() + Size * 1000;
                cromo = variant.getContig();
                array.add(num);
                num = 0;

            }
        }
        variante = "" + array.size();

    }
}
