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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Salma Alzahrani
 * This class is the core  part of PAREnet
 */
public class ParenetMain {

    /**
     * @param args the command line arguments
     */
    
    private final HashMap<String,Transcript> targetedTranscripts = new HashMap();
    private final PAREnet_NetworkConstructionInput input;
    
    
    public static void main(String[] args) {
        // TODO code application logic here
    }

    public ParenetMain() {
        System.out.println("----------------------------------");
        input = PAREnet_NetworkConstructionInput.getInstance();
    }

    
    public void process(){
         
        File ps2_out_file;

        if (input.getPAREsnip2OutputFile()!= null) {
            ps2_out_file = input.getPAREsnip2OutputFile(); 
            System.out.println("--- Skipping PAREsnip2 analysis ---");
        } else{
            ps2_out_file = new File(input.getOutputDir().getAbsolutePath() + File.separator + input.getSmallRNAFile().getName().split("\\.")[0] + "_"  + input.getDegradomeFile().getName().split("\\.")[0] + ".csv");
            // Start PAREsnip2 analysis
            System.out.println("----------------------------------");
            System.out.println("Start PAREsnip2 analysis.");
            System.out.println("----------------------------------");
            System.out.println();
            
            performPAREsnip2Analysis();

            System.out.println();
            System.out.println("----------------------------------");
            System.out.println("Finished PAREsnip2 analysis.");
            System.out.println("----------------------------------");
        }
        
        // 1. parse parsnip2 results
        // 2. find conservered (optional)
        // 2. get interactions
        // 3. set miRNAs annotation
        // 4. remove other RNAs
        // 5. Set source genes
        // 6. set target validation
        
        File functional_srnas_fasta = new File(ps2_out_file.getParent() + File.separator + ps2_out_file.getName().substring(0, ps2_out_file.getName().lastIndexOf(".")) + "_functional_srnas.fa");
        functional_srnas_fasta.deleteOnExit();
        
        Map<String, Set<Interaction>> interactions_srna_mrna = new HashMap();
        Map<String, Set<String>> interactions_mrna_srna = new HashMap();
        getPAREsnip2Results(interactions_srna_mrna, interactions_mrna_srna, ps2_out_file, true); // false: if test file (not ps2 file)
        
        Set<String> seqs = interactions_srna_mrna.keySet();
        writeToFasta(seqs,functional_srnas_fasta);
        
        HashMap<String,String> functional_mirnas = new HashMap();
        if(input.getMiRBaseFile()!=null){
            functional_mirnas = setMirnaAnnotation(interactions_srna_mrna, functional_srnas_fasta);
        }
        
        // Remove other RNAs (t-r-RNA) and their targeted transcripts
        //////////////////////////////////////////////////////////////   
        // Comment out if a correct directory of Rfam FASTA file is provided
//            Set<String> discarded_seqs = removeOtherRNAs(functional_srnas_fasta);
//            discarded_seqs.removeAll(functional_mirnas.keySet()); // to keep known miRNAs that are in Rfam
//
//            for(String seq : discarded_seqs){
//                for(Interaction inter : interactions_srna_mrna.get(seq)){
//                    this.targetedTranscripts.remove(inter.getmRNA().getGene_id());
//                    interactions_mrna_srna.remove(inter.getmRNA().getGene_id());
//                }
//            }
//            interactions_srna_mrna.keySet().removeAll(discarded_seqs);
        //////////////////////////////////////////////////////////////
        
        if(input.getMiRTarBaseFile() != null){
            setValidInteractions(interactions_srna_mrna, functional_mirnas);
        }
        
        writeToFasta(seqs,functional_srnas_fasta);
        
        setSourceGene(interactions_srna_mrna, interactions_mrna_srna,functional_srnas_fasta);
        
        File outFile = new File(input.getOutputDir() + File.separator + "PAREnet_CysIn.csv");
        writeResultsFile(interactions_srna_mrna, outFile);
        
        System.out.println();
        System.out.println("PAREnet network construction process is complete!");
        System.out.println("You can find the results in: " + outFile.getAbsolutePath());
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
    
    private  HashMap<String,String> setMirnaAnnotation(Map<String, Set<Interaction>> interactions, File srnasFile){
        System.out.println("Annotating valid miRNAs ...");
        File miRNAs = input.getMiRBaseFile();
        HashMap<String,String> functional_mirnas = new HashMap();

        File miRNAAlignment = patmanMapping(srnasFile, miRNAs, srnasFile.getParentFile().getAbsolutePath(), 0, 2);
        miRNAAlignment.deleteOnExit();
        HashSet<String> miRNAAlignedReads = getPatmanLines(miRNAAlignment);

        for(String line: miRNAAlignedReads){
            
            String splits[] = line.split("\t");
            String s = splits[1].split("\\(")[0].trim();

            String miRID = s;
            miRID = splits[0].split(" ")[0].substring(4).trim(); 

            functional_mirnas.put(s,miRID);
            for(Interaction i : interactions.get(s)){
                i.setMirnaAnnotation(miRID);
            }
        }
        return functional_mirnas;
    }
    
    private  void setSourceGene(Map<String, Set<Interaction>> interactions, Map<String, Set<String>> interactions_mrna_srna, File srnasFile){
        System.out.println("Setting source genes ...");
        File transcriptome = input.getTranscript();
        
        File sRNAAlignment = patmanMapping(srnasFile, transcriptome, srnasFile.getParentFile().getAbsolutePath(), 0, 0);
        sRNAAlignment.deleteOnExit();
        
        HashSet<String> srnaAlignedReads = getPatmanLines(sRNAAlignment);
        
        for(String line: srnaAlignedReads){
            String[] spl = line.split("\t");
            String[] splits = spl[0].split("\\|");
            
            String srna = spl[1].split("\\(")[0].trim();
            String gene_id = splits[0].split(" ")[0].trim();
            
            // Skip this transcript if it is not predicted in PAREsnip2 analysis
            if(!targetedTranscripts.containsKey(gene_id)){
                continue;
            }
            
            String mirid = srna;
            boolean ismirna = false;
            // To retrieve miRNA annotation
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
            // Update the targeted transcripts to be annotated as source
            for(String s : interactions_mrna_srna.get(gene_id)){
                for(Interaction i : interactions.get(s)){
                    if(i.getmRNA().getGene_id().contains(gene_id)){
                        i.getmRNA().setSource();
                    }
                }
            }
        }
    }
    
    
    private void setValidInteractions(Map<String, Set<Interaction>> interactions, HashMap<String,String> functional_mirnas ){
        System.out.println("Annotating experimentally validated targets ...");
        Map<String, List<String>> validInteractions = new HashMap();

        File validFile = input.getMiRTarBaseFile();
        readMitarbaseValidTargets(validFile, validInteractions);

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
                String mirid = splits[1].replace("ath-", "").trim();
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
  

    private static Set<String> removeOtherRNAs(File srnasFile) {
        
        System.out.println("Removing other non coding RNAs ...");
        File tRNA, rRNA;
        // Update directories
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
*/
        return patmanOutput;
    }
    
    // If PAREnip2 output file was provided
    private HashSet<String> getFunctionalsRNAs(File ps2_out_file){
        HashSet<String> functional_srnas = new HashSet();
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(ps2_out_file));
            String line;
            String srna;
            //header
            br.readLine();
            while ((line = br.readLine()) != null) {               
                srna = line.split("\"")[1].substring(1);
                functional_srnas.add(srna);
            }
        } catch (IOException ex) {
        }
        
        return functional_srnas;
    }
    
    // Paresing PAREsnip2 results
    private void getPAREsnip2Results(Map<String, Set<Interaction>> interactions_srna_mrna,Map<String, Set<String>> interactions_mrna_srna, File ps2_out_file, boolean isPS2File) {
        System.out.println("Parsing PAREsnip2 output table ...");

        try {
            BufferedReader br = new BufferedReader(new FileReader(ps2_out_file));
            String line;
            String srna;
            String gene_id;
            String gene_name; // gene symbol
            String gene_description;
            //skip header
            br.readLine();
            while ((line = br.readLine()) != null) {
                String[] fields1 = line.split("\"");
                String[] fields2 = fields1[4].split(",");
               
                srna = fields1[1].substring(1);
                int category = Integer.valueOf(fields2[1]);
                if(category>2){
                    br.readLine();
                    br.readLine();
                    continue;
                }
                int cleavePos = Integer.valueOf(fields2[2]);
                gene_description = fields1[3].replace(",", ";").trim();
                String fields3[] = gene_description.split("\\|");
                gene_id = fields3[0].split(" ")[0];
                gene_name = gene_id.split("\\.")[0].split("-")[0];     // set the gene name to gene id in case there's no defined symbol           

                // new transcript
                Transcript mRNA = new Transcript(gene_id, gene_name, gene_description);
                // check if new transcrpt already exist
                if(!targetedTranscripts.containsKey(gene_id)){
                    // Add new mRNA target
                    targetedTranscripts.put(gene_id, mRNA);
                    // extract the gene name/symbol from the gene description in the transcriptome fasta file
                    if(fields3.length>2){ // if symbol found
                        if(fields3[1].trim().split(" ").length > 1){
                        gene_name = fields3[1].split(":")[1].trim().replace(",", ";");
                        mRNA.setGeneName(gene_name);
                        }
                    } else if(fields3.length>1){ // if no symbol found
                        // set the gene name with gene id 
                        mRNA.setGeneName(gene_name);
                    }
                }else {
                    gene_name = mRNA.getGene_name();
                }
                 
                Interaction i = new Interaction(srna,category, cleavePos, mRNA);
                // if sRNA is parsed for the first time
                if (!interactions_srna_mrna.containsKey(srna)) {
                    interactions_srna_mrna.put(srna, new HashSet());
                }
                // if this is the first occurance of interaction
                if(!interactions_srna_mrna.get(srna).contains(i)){
                    interactions_srna_mrna.get(srna).add(i);
                }
                // if this is the first occurance of mRNA
                if (!interactions_mrna_srna.containsKey(gene_id)) {
                    interactions_mrna_srna.put(gene_id, new HashSet());
                }
                if(!interactions_mrna_srna.get(gene_id).contains(srna)){
                    interactions_mrna_srna.get(gene_id).add(srna);
                }
                // skip two lines to get to the next interaction
                if(isPS2File){
                    br.readLine();
                    br.readLine();
                }
            }
        } catch (IOException | NumberFormatException ex) {
        }
    }
    
  
    static Set<String> getSeqsFasta(File seqsFile) {

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
    
     private static void writeResultsFile( Map<String, Set<Interaction>> interactions, File toWrite) {
         System.out.println("Writing results to output file ...");
         // Results header
         String header = ("source_node,sRNA seq,is miRNA,node_label,node_color,target_node,Gene description,node_label,node_color,Is source gene(node),Is source gene(edge),is validated target,Category,Cleavage position");

         try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(toWrite));

            bw.write(header + "\n"); // new line
            // for each sRNA 
            for (String seq : interactions.keySet()) {
                // for each interaction for that sRNA
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
     
     private void performPAREsnip2Analysis() {
/*
        // Initiliasing and setting PAREsnip2 configuration 
        Paresnip2Configuration ps2Config = Paresnip2Configuration.getInstance();
        RuleSet ps2Ruleset = RuleSet.getRuleSet();
        // Initiliasing PAREsnip2 input
        Paresnip2InputFiles ps2_input = Paresnip2InputFiles.getInstance();
            
        if (input.getPAREsnip2Config() == null) {
            ps2Config.setDefaultNATsiParameters();
        } else {
            try {
                ps2Config.loadConfig(input.getPAREsnip2Config());
            } catch (Exception ex) {
                Logger.getLogger(ParenetMain.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        // if Targetting rule set, the default target (Allen rules) is selected
        if (input.getPAREsnip2TargetingRules() == null) {
            ps2Ruleset.setDefaultAllen();
        } else {
            try {
                ps2Ruleset.loadRules(input.getPAREsnip2TargetingRules());
            } catch (IOException | NumberFormatException ex) {
                Logger.getLogger(ParenetMain.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        // Adding input files
        ps2_input.addSmallRNAReplicate(input.getSmallRNAFile());
        ps2_input.addTranscriptome(input.getTranscript());
        ps2_input.addDegradomeReplicate(input.getDegradomeFile());
        ps2_input.setOuputDirectory(input.getOutputDir());
            
        // check if GFF file is provided
        if (ps2Config.isUsingGenomeGFF()) {
            if(input.getGFF_file()==null){
                System.out.println("ERROR - You configuration is set to use a genome + corresponding gff3 file but the -gff3 path-to-file parameter is missing");
                System.exit(1);
            }
            ps2_input.addGFF(input.getGFF_file());
            ps2_input.addGenome(input.getTranscript());
        }
        // check if genome file is provided
        if (Paresnip2Configuration.getInstance().isAlignSmallRNAsToGenome()) {
            if(input.getGenomeFile()==null){
            System.out.println("ERROR - You configuration is set to use a genome + corresponding gff3 file but the -gff3 path-to-file parameter is missing");
            System.exit(1);
        }
            ps2_input.addGenome(input.getGenomeFile());
        }
        
        Engine ps2_engine = new Engine();
*/
    }
    
    
}
