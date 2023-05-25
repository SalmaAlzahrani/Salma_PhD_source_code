/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package parenet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Salma Alzahrani
 * 
 * This class is for testing purposes
 */
public class main {

    private HashMap<String,Transcript> targetedTranscripts = new HashMap();

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
    }
    
     public void startProcess(){
/*
//        File ps2_out_rep1 = new File("/Users/salmayz/Documents/networks_project_start/reads/A/carrington_rules/no_muliplier/cat_0-3/conserved/WTA_Col_A.csv"); // Terst file
//        File ps2_out_rep2 = new File("/Users/salmayz/Documents/networks_project_start/reads/A/carrington_rules/no_muliplier/cat_0-3/conserved/WTA_Col_B.csv"); // Terst file
//        File ps2_out_rep3 = new File("/Users/salmayz/Documents/networks_project_start/reads/A/carrington_rules/no_muliplier/cat_0-3/conserved/WTA_Col_C.csv"); 
        
        File ps2_out_rep1 = new File("/Users/salmayz/Documents/networks_project_start/reads/A/carrington_rules/no_muliplier/cat_0-2/WTA_Col_A_ACAGTG_L004_R1.csv");
        File ps2_out_rep2 = new File("/Users/salmayz/Documents/networks_project_start/reads/B/carrington_rules/no_muliplier/cat_0-2/WTB_Col_B_ACTTGA_L004_R1.csv");
        File ps2_out_rep3 = new File("/Users/salmayz/Documents/networks_project_start/reads/C/carrington_rules/no_muliplier/cat_0-2/WTC_Col_C_GTAGAG_L004_R1.csv");

//        File ps2_out_rep1 = new File("/Users/salmayz/Documents/networks_project_start/reads/A/allen_rules/cat0-3_noMul/conserved/WTA_Col_A.csv");
//        File ps2_out_rep2 = new File("/Users/salmayz/Documents/networks_project_start/reads/A/allen_rules/cat0-3_noMul/conserved/WTB_Col_B.csv");
//        File ps2_out_rep3 = new File("/Users/salmayz/Documents/networks_project_start/reads/A/allen_rules/cat0-3_noMul/conserved/WTC_Col_C.csv");


        File functional_srnas_fasta = new File(ps2_out_rep1.getParent() + File.separator + ps2_out_rep1.getName().substring(0, ps2_out_rep1.getName().lastIndexOf(".")) + "_functional_srnas_test.fa");
        functional_srnas_fasta.deleteOnExit();
        

        // Conservation 2 out of 3 reps
//        File conservedFile = new File(ps2_out_rep1.getParent() + File.separator + ps2_out_rep1.getName().substring(0, ps2_out_rep1.getName().lastIndexOf(".")) +  "_conserved_2_out_3_2022-711.csv");
//        findOverlapBetweenTargetPrediction( ps2_out_rep1,true, ps2_out_rep2, true, conservedFile, true);
//        findOverlapBetweenTargetPrediction( ps2_out_rep2, true,  ps2_out_rep3, true, conservedFile, false);
//        findOverlapBetweenTargetPrediction( ps2_out_rep3, true, ps2_out_rep1, true, conservedFile, false);
        
        // Conservation 3 out of 3 reps
        File conservedFile = new File(ps2_out_rep1.getParent() + File.separator + ps2_out_rep1.getName().substring(0, ps2_out_rep1.getName().lastIndexOf(".")) + "_C_conserved_3_out_3.csv");
        conservedFile.deleteOnExit();
        findOverlapBetweenTargetPrediction( ps2_out_rep1, true, ps2_out_rep2, true, conservedFile, true);
        findOverlapBetweenTargetPrediction( ps2_out_rep3, true, conservedFile, false, conservedFile, true);
        
        
//         Map<String, Set<Interaction>> interactions = getPAREsnip2Results(ps2_out_rep1,true);
        Map<String, Set<Interaction>> interactions = getPAREsnip2Results( conservedFile, false);
      
        
        Set<String> seqs = interactions.keySet();
        writeToFasta(seqs,functional_srnas_fasta);
        HashMap<String,String> functional_mirnas = setMirnaAnnotation(interactions, functional_srnas_fasta);
        
        
        // Remove other RNAs (t-r-RNA) and their targeted transcripts .. REMOVE OTHER SRNAS THEN SET SOURCE GENES
        Set<String> discarded_seqs = removeOtherRNAs(functional_srnas_fasta);
        discarded_seqs.removeAll(functional_mirnas.keySet()); // to keep known miRNAs that are in Rfam
        interactions.keySet().removeAll(discarded_seqs);
        for(String seq : discarded_seqs){
            for(Interaction inter : interactions.get(seq)){
                this.targetedTranscripts.remove(inter.getGene_ID());
            }
        }

        setValidInteractions(interactions, functional_mirnas);
        
        writeToFasta(seqs,functional_srnas_fasta);
        setSourceGene(interactions,functional_srnas_fasta);
        // HashSet<String> sourceGenes = getSourceGenes(interactions, conservedFile, functional_srnas_fasta);
        
        // File outFile = new File(ps2_out_rep1.getParent()+ File.separator + ps2_out_rep1.getName().substring(0, ps2_out_rep1.getName().lastIndexOf(".")) + "fileid.CysIn.csv");
        File outFile = new File(conservedFile.getParent() + File.separator + conservedFile.getName().substring(0, conservedFile.getName().lastIndexOf(".")) + "_withSource_20220808_test.CysIn.csv");

        writeResultsFile(interactions, outFile);
*/
    }
     

     private void getPAREsnip2Results(Map<String, Set<Interaction>> interactions_srna_mrna,Map<String, Set<String>> interactions_mrna_srna, File ps2_out_file, boolean isPS2File) {
        System.out.println("Parsing PAREsnip2 output table ...");

        try {
            BufferedReader br = new BufferedReader(new FileReader(ps2_out_file));
            String line;
            String srna;
            String gene_id;
            String gene_name; // Symbol
            String gene_description;
            //boolean isSourceGene;
            //header
            br.readLine();
            while ((line = br.readLine()) != null) {
                //isSourceGene = false;

                String[] fields1 = line.split("\"");
                String[] fields2 = fields1[4].split(",");
               
                srna = fields1[1].substring(1);
                int category = Integer.valueOf(fields2[1]);
                if(isPS2File){
                    br.readLine();
                    br.readLine();
                }
                
                int cleavePos = Integer.valueOf(fields2[2]);
                
                gene_description = fields1[3].replace(",", ";").trim();
                String fields3[] = gene_description.split("\\|");
                gene_id = fields3[0].split(" ")[0];
                gene_name = gene_id.split("\\.")[0];                

                Transcript mRNA = new Transcript(gene_id, gene_name, gene_description);
                if(!targetedTranscripts.containsKey(gene_id)){
                    targetedTranscripts.put(gene_id, mRNA);
                }
                if(fields3[1].trim().split(" ").length > 1){
                    gene_name = fields3[1].split(":")[1].trim().replace(",", ";");
                    mRNA.setGeneName(gene_name);
                }
                
                Interaction i = new Interaction(srna,category, cleavePos, mRNA);

                if (!interactions_srna_mrna.containsKey(srna)) {
                    interactions_srna_mrna.put(srna, new HashSet());
                }
                if(!interactions_srna_mrna.get(srna).contains(i)){
                    interactions_srna_mrna.get(srna).add(i);
                }
                
                if (!interactions_mrna_srna.containsKey(gene_id)) {
                    interactions_mrna_srna.put(gene_id, new HashSet());
                }
                if(!interactions_mrna_srna.get(gene_id).contains(srna)){
                    interactions_mrna_srna.get(gene_id).add(srna);
                }
            }
        } catch (IOException | NumberFormatException ex) {
        }

    }
    
    
    private static Map<String, List<String>> getPAREsnip2Results1(File results, boolean isPS2results) {
        Map<String, List<String>> targets = new HashMap();

        try {
            BufferedReader br = new BufferedReader(new FileReader(results));
            String line;
            //header
            br.readLine();
            while ((line = br.readLine()) != null) {
                String splits[] = line.split(",");
                String seq = splits[1].replace("\"", "").substring(1);
                String gene = splits[2].replace("\"", "").split(" ")[0];

                if (!targets.containsKey(seq)) {
                    targets.put(seq, new ArrayList());
                }

                if (!targets.get(seq).contains(gene)) {
                    targets.get(seq).add(gene);
                }
                if(isPS2results){
                    br.readLine();
                    br.readLine();
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return targets;
    }
     
     // Function needed when adding conservation feature to the tool
    private static void findOverlapBetweenTargetPrediction( File rep1, boolean isRep1PS2, File rep2, boolean isRep2PS2, File conservedFile, boolean printHeader) {

        String header = "";
        Map<String, List<String>> rep1Targets = getPAREsnip2Results1(rep1,isRep1PS2);
        Map<String, List<String>> rep2Targets = getPAREsnip2Results1(rep2,isRep2PS2);
        
        System.out.println( rep1.getName() + " : " + rep1Targets.size());
        System.out.println(rep1.getName() + " : " + rep2Targets.size());
        
        Map<String, List<String>> conserved = new HashMap();

        for (String seq : rep1Targets.keySet()) {
            if (rep2Targets.containsKey(seq)) {
                for (String target : rep1Targets.get(seq)) {
                    if (rep2Targets.get(seq).contains(target)) {
                        if (!conserved.containsKey(seq)) {
                            conserved.put(seq, new ArrayList<>());
                        }

                        conserved.get(seq).add(target);
                    }
                }
            }
        }

        System.out.println("Conserved interactions: " + conserved.size());
        
        Set<String> toWrite = new HashSet();

        try {
            BufferedReader br = new BufferedReader(new FileReader(rep1));
            String line;
            header = br.readLine();
            while ((line = br.readLine()) != null) {
                line += br.readLine();
                line += br.readLine();

                String splits[] = line.split(",");
                String sRNA = splits[1].replace("\"", "").substring(1);
                String gene = splits[2].replace("\"", "").split(" ")[0];

                if (conserved.containsKey(sRNA) && conserved.get(sRNA).contains(gene)) {
                    toWrite.add(line);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
        System.out.println("to write interactions: " + toWrite.size());
        
        //Go through the results file and spit out the conserved
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(conservedFile, !printHeader));
            if(printHeader)
                bw.write(header + "\n");
            for (String result : toWrite) {
                bw.write(result + "\n");
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println();
    }
     
    
    private  void setSourceGene(Map<String, Set<Interaction>> interactions, Map<String, Set<String>> interactions_mrna_srna, File srnasFile){
        System.out.println("Setting source genes ...");
        File transcriptome = new File("");
        
        File sRNAAlignment = patmanMapping(srnasFile, transcriptome, srnasFile.getParentFile().getAbsolutePath(), 0, 0);
        sRNAAlignment.deleteOnExit();
        
        HashSet<String> srnaAlignedReads = getPatmanLines(sRNAAlignment);
        
        for(String line: srnaAlignedReads){
            String[] spl = line.split("\t");
            String[] splits = spl[0].split("\\|");
            
            String srna = spl[1].split("\\(")[0].trim();
            String gene_id = splits[0].trim();
            
            if(!targetedTranscripts.containsKey(gene_id)){
                continue;
            }
            
            String mirid = srna;
            boolean ismirna = false;
            // To get mirna annotation
            for(Interaction i : interactions.get(srna)){
                mirid = i.getMiRid();
                ismirna = i.isIsKnownMirna();
                break;
        }
            Transcript mRNA = targetedTranscripts.get(gene_id);
            mRNA.setSource();

            Interaction new_i = new Interaction(srna, mirid, ismirna, mRNA);

            if(!interactions.get(srna).contains(new_i)){
                interactions.get(srna).add(new_i);
            }
            
            for(String s : interactions_mrna_srna.get(gene_id)){
                for(Interaction i : interactions.get(s)){
                    if(i.getmRNA().getGene_id().contains(gene_id)){
                        i.getmRNA().setSource();
                    }
                }
            }
        }
    }
    
     private  HashMap<String,String> setMirnaAnnotation(Map<String, Set<Interaction>> interactions, File srnasFile){
        File miRNAs = new File("/Users/salmayz/Documents/networks_project_start/mature.fa");
        HashMap<String,String> functional_mirnas = new HashMap();

        File miRNAAlignment = patmanMapping(srnasFile, miRNAs, srnasFile.getParentFile().getAbsolutePath(), 0, 2);
        miRNAAlignment.deleteOnExit();
        HashSet<String> miRNAAlignedReads = getPatmanLines(miRNAAlignment);

        for(String line: miRNAAlignedReads){
            
            String splits[] = line.split("\t");
            String s = splits[1].split("\\(")[0].trim();

            String miRID = s;
            miRID = splits[0].split(" ")[0].replace("ath-", "").trim();

            functional_mirnas.put(s,miRID);
            for(Interaction i : interactions.get(s)){
                i.setMirnaAnnotation(miRID);
            }
        }
        
        return functional_mirnas;
    }
    
    private static Set<String> removeOtherRNAs(File srnasFile) {
        File tRNA, rRNA;

        tRNA = new File("/Users/salmayz/WorkbenchFull/data/t_and_r_RNAs.fa");
        rRNA = new File("/Users/salmayz/WorkbenchFull/data/Rfam.fasta");

        Set<String> seqs = new HashSet();

        
        File tRNAAlignment = patmanMapping(srnasFile, tRNA, srnasFile.getParentFile().getAbsolutePath(), 0, 0);
        tRNAAlignment.deleteOnExit();
        File rRNAAlignment = patmanMapping(srnasFile, rRNA, srnasFile.getParentFile().getAbsolutePath(), 0, 0);
        rRNAAlignment.deleteOnExit();

        
        HashSet<String> tRNAAlignedReads = getPatmanSeqs(tRNAAlignment);
        HashSet<String> rRNAAlignedReads = getPatmanSeqs(rRNAAlignment);

        seqs.addAll(tRNAAlignedReads);
        seqs.addAll(rRNAAlignedReads);

        return seqs;
    }
    
    private static HashSet<String> getPatmanSeqs(File alignmentFile) {
        HashSet<String> alignedReads = new HashSet();

        try {
            BufferedReader br = new BufferedReader(new FileReader(alignmentFile));
            String line;
            while ((line = br.readLine()) != null) {
                String splits[] = line.split("\t");
                String seq = splits[1].split("\\(")[0].trim();
                alignedReads.add(seq);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return alignedReads;
    }
    
    private void setValidInteractions(Map<String, Set<Interaction>> interactions, HashMap<String,String> functional_mirnas ){
        Map<String, List<String>> validInteractions = new HashMap();

        File validFile = new File("/Users/salmayz/Documents/networks_project_start/ath_MTI.csv");
        readMitarbaseValidTargets(validFile, validInteractions);
        validFile = new File("/Users/salmayz/Documents/networks_project_start/validated_targets_psRNATarget.csv");
        readPsrnaValidTargets(validFile, validInteractions);
        validFile = new File("/Users/salmayz/Documents/networks_project_start/reads/validated_targets.csv");
        readOtherValidTargets(validFile, validInteractions, interactions);
        
        for(String srna: functional_mirnas.keySet()){

            String mirid = functional_mirnas.get(srna); 
            if(validInteractions.containsKey(mirid)){
                for(Interaction i : interactions.get(srna)){
                    for(String validGene : validInteractions.get(mirid)){
                        String geneName = i.getGene_name();
                        String geneId = i.getGene_ID();
                        if(geneName.contains(validGene) || geneId.contains(validGene)){
                            i.setValidated();
                        }
                    }
                }
            }
        }
    }
    
    private void readMitarbaseValidTargets(File validFile , Map<String, List<String>> validInteractions){
        try {
            BufferedReader br = new BufferedReader(new FileReader(validFile));
            br.readLine(); // header
            String line;
            while ((line = br.readLine()) != null) {
                String splits[] = line.split(",");
                // MIRT002091,ath-miR398c-3p,Arabidopsis thaliana,CSD2,817365,Arabidopsis thaliana,"5""RACE//Northern blot",Functional MTI (Weak),20400846
                String mirid = splits[1].replace("ath-", "").trim(); //split("-")[1].trim();;
                String gene = splits[3].trim();
                if(!validInteractions.containsKey(mirid)){
                    validInteractions.put(mirid, new ArrayList<>());
                }
                if (!validInteractions.get(mirid).contains(gene)) {
                    validInteractions.get(mirid).add(gene);
                }
                
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    
    private void readPsrnaValidTargets(File validFile , Map<String, List<String>> validInteractions){
        try {
            BufferedReader br = new BufferedReader(new FileReader(validFile));
            br.readLine();
            br.readLine(); // header
            String line;
            while ((line = br.readLine()) != null) {
                String splits[] = line.split(",");
                String mirid = splits[0].replace("ath-", "").trim(); //split("-")[1].trim();;
                String gene = splits[2].trim();
                if(!validInteractions.containsKey(mirid)){
                    validInteractions.put(mirid, new ArrayList<>());
                }
                if (!validInteractions.get(mirid).contains(gene)) {
                    validInteractions.get(mirid).add(gene);
                }
                
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private void readOtherValidTargets(File validFile , Map<String, List<String>> validInteractions, Map<String, Set<Interaction>> interactions){
         try {
            BufferedReader br = new BufferedReader(new FileReader(validFile));
            br.readLine(); // header
            String line;
            while ((line = br.readLine()) != null) {
                String splits[] = line.split(",");
                String mirid_splits[] = splits[0].substring(1).split(">");                
                //String mirid = splits[1].split("-")[1].trim();
                String gene = splits[2].trim();
                for(String mirid : mirid_splits){
                    mirid = mirid.replace("ath-", "").trim(); //split("-")[1].trim();
                    if(!validInteractions.containsKey(mirid)){
                        validInteractions.put(mirid, new ArrayList<>());
                    }
                    if (!validInteractions.get(mirid).contains(gene)) {
                        validInteractions.get(mirid).add(gene);
                    }
                }
                
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
//    private void readSmallRNAs() {
//
//        File file = new File("/Users/salmayz/Transcriptome/TAIR10_cdna_20110103_representative_gene_model_updated_filtered_MT_C.fa");
//        File tutorial = new File("/Users/salmayz/WorkbenchFull/PAREnet_test_files/transcript.fa");
//
//        File out = new File("/Users/salmayz/WorkbenchFull/PAREnet_test_files/transcript_update.fa");
//
//        HashSet<String> smallRNAs = new HashSet();
//        String line;
//
//        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
//
//            while ((line = br.readLine()) != null) {
//                    if (line.startsWith(">")) {
//                        
//                        smallRNAs.add(line.trim());
//                    }
//                
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//
//        try {
//            BufferedWriter bw = new BufferedWriter(new FileWriter(out));
//            try (BufferedReader br = new BufferedReader(new FileReader(tutorial))) {
//
//                while ((line = br.readLine()) != null) {
//                        if (line.startsWith(">")) {
//                            for(String s: smallRNAs){
//                                if (s.contains(line)) {
//                                    bw.write(s + "\n");
//                                    bw.write(br.readLine() + "\n");
//                                }
//                            }
//                        }
//                }
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//
//            bw.flush();
//            bw.close();
//
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//        
//    }
  
    private static File patmanMapping(File sRNAs, File longReads, String outputDirectory, int gaps, int mismatches) {

        String alignedFile = outputDirectory + File.separator + sRNAs.getName() + "_" + longReads.getName() + ".patman";

        File patmanOutput = new File(alignedFile);
/*
        PatmanParams newP_Params = new PatmanParams();
        newP_Params.setMaxGaps(gaps);
        newP_Params.setMaxMismatches(mismatches);
        newP_Params.setPreProcess(false);
        newP_Params.setPostProcess(false);
        newP_Params.setMakeNR(true);
        newP_Params.setPositiveStrandOnly(false);

        File tmpDir = new File(outputDirectory + File.separator + "/tmp");
        tmpDir.deleteOnExit();
        try {
            tmpDir.mkdir();
            PatmanRunner runner = new PatmanRunner(sRNAs, longReads,
                    patmanOutput, tmpDir, newP_Params);
            runner.setUsingDatabase(false);

            Thread myThread = new Thread(runner);

            myThread.start();
            myThread.join();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        tmpDir.deleteOnExit();
*/
        return patmanOutput;
    }
    
    private static HashSet<String> getPatmanLines(File alignmentFile) {
        HashSet<String> line_seq =  new HashSet();
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(alignmentFile));
            String line;
            while ((line = br.readLine()) != null) {
                line_seq.add(line);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
        return line_seq;
    }
    
    
    private static Set<String> getSeqsFasta(File seqsFile) {

        HashSet<String> seqs = new HashSet();
        try {
            BufferedReader br = new BufferedReader(new FileReader(seqsFile));
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith(">")) {
                    seqs.add(line);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return seqs;
    }
    
    private void writeToFasta(Set<String> seqs, File file){
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));

            for (String seq : seqs) {
                bw.write(">" + seq + "\n");
                bw.write(seq + "\n");
            }

            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
     private static void writeResultsFile( Map<String, Set<Interaction>> interactions, File toWrite) {
        String header = ("miRid(Source node),sRNA seq,is miRNA,node_color,node_color,Gene name(target node),Gene ID,Gene description,Is source gene(node),Is source gene(edge),is validated target,Category,Cleavage position");

         try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(toWrite));

            bw.write(header + "\n");
            for (String seq : interactions.keySet()) {
                for(Interaction i : interactions.get(seq)){
                    bw.write(i.toString()+ "\n");
                }
            }

            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }
     
       /*
    private HashSet<String> getSourceGenes(Map<String, Set<Interaction>> interactions, File ps2File, File srnasFile){
        Map<String, List<String>> targets = getPAREsnip2Results1(ps2File,false);
        writeToFasta(targets.keySet(), srnasFile);
        
        Set<String> targetedGenes = new HashSet();
        
        for(String s: targets.keySet()){
            for(String g : targets.get(s)){
                targetedGenes.add(g);
            }
        }
        
        File transcriptome = new File("/Users/salmayz/Documents/networks_project_start/Transcriptome/TAIR10_cdna_20110103_representative_gene_model_updated_filtered_MT_C.fa");
        HashSet<String> sourceGenes = new HashSet();
        
        File sRNAAlignment = patmanMapping(srnasFile, transcriptome, srnasFile.getParentFile().getAbsolutePath(), 0, 0);
        sRNAAlignment.deleteOnExit();
        
        HashSet<String> srnaAlignedReads = getPatmanLines(sRNAAlignment);
        
        for(String line: srnaAlignedReads){
            String[] spl = line.split("\t");
            String[] splits = spl[0].split("\\|");
         
            String srna = spl[1].split("\\(")[0].trim();
            String gene_id = splits[0].trim();

            if(targetedGenes.contains(gene_id)){
                String gene_name = gene_id;
                String gene_description = spl[0].replace(",", ";");
            
                sourceGenes.add(gene_id);
            
                String mirid = srna;
                boolean ismirna = false;

                Interaction new_i = new Interaction(srna, mirid, gene_id, gene_name, gene_description, ismirna);

                if(splits[1].trim().split(" ").length > 1){
                    gene_name = splits[1].split(":")[1].trim().replace(",", ";");
                    new_i.setGeneName(gene_name);
                }
                
                if (!interactions.containsKey(srna)) {
                    interactions.put(srna, new HashSet());
                }
                if(!interactions.get(srna).contains(new_i)){
                    interactions.get(srna).add(new_i);
                }
            }
        }
        
        return sourceGenes;
    }
    */  
     
     
      /*
     private static void readAndWriteTestFile(){
        File inFile = new File ("/Users/salmayz/Documents/networks_project_start/reads/output_file_20220615_2.csv");
        File outFileCys = new File(inFile.getParent() + File.separator + inFile.getName().substring(0, inFile.getName().lastIndexOf(".")) + ".removed_t_r_RNA.CysIn.csv");
        File functional_srnas_fasta = new File(inFile.getParent() + File.separator + inFile.getName().substring(0, inFile.getName().lastIndexOf(".")) + "_functional_srnas.fa");
        functional_srnas_fasta.deleteOnExit();
        Map<String, Set<Interaction>> interactions = new HashMap<>();
        
        // parse inFile
        try {
            BufferedReader br = new BufferedReader(new FileReader(inFile));
            String line;
            String srna;
            String gene_id;
            String gene_name; // Symbol
            String gene_description;

            br.readLine();  // header
            while ((line = br.readLine()) != null) {

                String[] fields = line.split(",");
               
                srna = fields[2].trim();
                int catA = Integer.valueOf(fields[3]);
                int catB = Integer.valueOf(fields[4]);
                int catC = Integer.valueOf(fields[5]);
                
                int cat = catA;
                if(catB < cat)
                    cat = catB;
                if(catC < cat)
                    cat = catC;
                int cleavePos = 0;
                boolean isSourceGene = false;
                gene_id = fields[1].trim();
                gene_name = fields[0].trim();
                gene_description = "";
                
                Interaction i = new Interaction(srna,  gene_id,  gene_name, gene_description, cat, cleavePos, isSourceGene);
                if(srna.contains("AAAATGGTGGAGAAAGAAGAG") || srna.contains("AAAAGGTTTTTCTTAGAGGGAC")){
                    System.out.println("");
                }
                if (!interactions.containsKey(srna)) {
                    interactions.put(srna, new HashSet());
                }
                if(!interactions.get(srna).contains(i)){
                    interactions.get(srna).add(i);
                }

            }
        } catch (IOException | NumberFormatException ex) {
        }
        
        writeToFasta(interactions.keySet(),functional_srnas_fasta);
        HashMap<String,String> functional_mirnas = setMirnaAnnotation(interactions, functional_srnas_fasta);
        // remove other rnas        
        Set<String> discarded_seqs = removeOtherRNAs(functional_srnas_fasta);
        discarded_seqs.removeAll(functional_mirnas.keySet()); // to keep known miRNAs that are in Rfam
        interactions.keySet().removeAll(discarded_seqs);
 
        setValidInteractions(interactions, functional_mirnas);
        writeResultsFile(interactions, outFileCys);
     }
     
*/
}
