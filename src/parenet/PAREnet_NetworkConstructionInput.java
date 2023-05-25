/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package parenet;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
//import javax.json.JsonObject;
//import uk.ac.uea.cmp.srnaworkbench.utils.JsonUtils;

/**
 *
 * @author Salma Alzahrani
 */
public class PAREnet_NetworkConstructionInput {
    
    private static PAREnet_NetworkConstructionInput input = null;
    File outputDir;
    File paresnip2OutFile;
    File smallRNAFile;
    File transcript;
    File degradomeFile;
    File GFF_file;
    File PAREsnip2Config;
    File PAREsnip2TargetingRules;
    File genomeFile;
    File mirtarbaseFile;
    File mirbaseFile;

    public static PAREnet_NetworkConstructionInput getInstance() {
        if (input == null) {
            input = new PAREnet_NetworkConstructionInput();
        }

        return input;
    }

    public File getOutputDir() {
        return outputDir;
    }

    public void setOutputDir(File outputDir) {

        if (!outputDir.exists()) {
            outputDir.mkdirs();
        }

        this.outputDir = outputDir;
    }

    public File getPAREsnip2OutputFile() {
        return paresnip2OutFile;
    }
    
    public File getSmallRNAFile() {
        return smallRNAFile;
    }

    public void setSmallRNAFile(File smallRNAFile) {
        this.smallRNAFile = smallRNAFile;
    }

    public File getTranscript() {
        return transcript;
    }

    public void setTranscript(File transcript) {
        this.transcript = transcript;
    }

    public File getDegradomeFile() {
        return degradomeFile;
    }

    public void setDegradomeFile(File degradomeFile) {
        this.degradomeFile = degradomeFile;
    }

    public File getGFF_file() {
        return GFF_file;
    }

    public void setGFF_file(File GFF_file) {
        this.GFF_file = GFF_file;
    }

    public File getPAREsnip2Config() {
        return PAREsnip2Config;
    }

    public void setPAREsnip2Config(File PAREsnip2Config) {
        this.PAREsnip2Config = PAREsnip2Config;
    }

    public File getPAREsnip2TargetingRules() {
        return PAREsnip2TargetingRules;
    }

    public void setPAREsnip2TargetingRules(File PAREsnip2TargetingRules) {
        this.PAREsnip2TargetingRules = PAREsnip2TargetingRules;
    }
    
    public File getMiRBaseFile() {
        return mirbaseFile;
    }
    
    public File getMiRTarBaseFile() {
        return mirtarbaseFile;
    }

    void reset() {
        input = null;
    }

    public File getGenomeFile() {
        return genomeFile;
    }

    public void setGenomeFile(File genomeFile) {
        this.genomeFile = genomeFile;
    }

    public void loadFromJson(File inputJson) throws IOException {
        /*
        JsonObject jsonObject = null;

        jsonObject = JsonUtils.parseJsonFile(inputJson);

        if (jsonObject.containsKey("paresnip2_output_file")) {
            Path paresnip2_output_file_path = Paths.get(jsonObject.getString("paresnip2_output_file"));
            if (Files.exists(paresnip2_output_file_path)) {
                paresnip2OutFile = paresnip2_output_file_path.toFile();
                System.out.println("PAREsnip2 output file is loaded.");
            } else {
                throw new IOException("PAREsnip2 output file: " + paresnip2_output_file_path + " not found.");
            } 
            
            if (jsonObject.containsKey("transcriptome_file")) {
                Path transcriptome_file_Path = Paths.get(jsonObject.getString("transcriptome_file"));
                if (Files.exists(transcriptome_file_Path)) {
                    transcript = transcriptome_file_Path.toFile();
                    System.out.println("A transcriptome file is loaded.");
                } else {
                    throw new IOException("transcriptome file: " + transcriptome_file_Path + " not found.");
                }
            } else {
                throw new IOException("A transcriptome file is required.");
            }
            
        } else if(jsonObject.containsKey("srna_file")){
            Path srna_file_Path = Paths.get(jsonObject.getString("srna_file"));
            if (Files.exists(srna_file_Path)) {
                smallRNAFile = srna_file_Path.toFile();

            } else {
                throw new IOException("sRNA file: " + srna_file_Path + " not found.");
            }
            System.out.println("A sRNA file is loaded.");       
            if (jsonObject.containsKey("transcriptome_file")) {
                Path transcriptome_file_Path = Paths.get(jsonObject.getString("transcriptome_file"));
                if (Files.exists(transcriptome_file_Path)) {
                    transcript = transcriptome_file_Path.toFile();
                    System.out.println("A transcriptome file is loaded.");

                } else {
                    throw new IOException("transcriptome file: " + transcriptome_file_Path + " not found.");
                }
                
            } else if ((jsonObject.containsKey("gff3_file") && jsonObject.containsKey("genome_file"))) {
                //gff
                Path gff_file_Path = Paths.get(jsonObject.getString("gff3_file"));
                if (Files.exists(gff_file_Path)) {
                    GFF_file = gff_file_Path.toFile();

                } else {
                    throw new IOException("GFF3 file: " + gff_file_Path + " not found.");
                }

                //genome
                Path genome_file_Path = Paths.get(jsonObject.getString("genome_file"));
                if (Files.exists(genome_file_Path)) {
                    genomeFile = genome_file_Path.toFile();

                } else {
                    throw new IOException("Genome file: " + genome_file_Path + " not found.");
                }

            } else {
                throw new IOException("You must provide a reference transcriptome in either FASTA or GFF3 format (with genome).");
            }

            if (jsonObject.containsKey("degradome_file")) {
                Path degradome_file_Path = Paths.get(jsonObject.getString("degradome_file"));
                if (Files.exists(degradome_file_Path)) {
                    degradomeFile = degradome_file_Path.toFile();
                    System.out.println("A degradome file is loaded.");
                } else {
                    throw new IOException("Degradome file: " + degradome_file_Path + " not found.");
                }
            }else{
                throw new IOException("A degradome file is required.");
            }
            //get PAREsnip2 config/target rules. If not present, use default.
            if (jsonObject.containsKey("paresnip2_config")) {
                Path paresnip2_config_file_Path = Paths.get(jsonObject.getString("paresnip2_config"));
                if (Files.exists(paresnip2_config_file_Path)) {
                    PAREsnip2Config = paresnip2_config_file_Path.toFile();

                } else {
                    throw new IOException("PAREsnip2 config file: " + paresnip2_config_file_Path + " not found.");
                }
            } else{
                System.out.println("PAREsnip2 config not set, using default configuration.");
            }

            if (jsonObject.containsKey("paresnip2_rules")) {
                Path paresnip2_rules_file_Path = Paths.get(jsonObject.getString("paresnip2_rules"));
                if (Files.exists(paresnip2_rules_file_Path)) {
                    PAREsnip2Config = paresnip2_rules_file_Path.toFile();

                } else {
                    throw new IOException("PAREsip2 rules file: " + paresnip2_rules_file_Path + " not found.");
                }
            } else {
                System.out.println("PAREsnip2 targeting criteria not set, using default configuration.");
            }
            
        }else {
            throw new IOException("You must provide a sRNA file or PAREsnip2 output file and a reference transcriptome");
        }
        
        
        if (jsonObject.containsKey("validated_targets_file")) {
            Path valid_targets_file_Path = Paths.get(jsonObject.getString("validated_targets_file"));
            if (Files.exists(valid_targets_file_Path)) {
                mirtarbaseFile = valid_targets_file_Path.toFile();
                System.out.println("Experimentally validated sRNA targets file (from miRTarBase) is loaded.");
            } else {
                throw new IOException("Experimentally validated sRNA targets file (from miRTarBase): " + valid_targets_file_Path + " not found.");
            }
            if (jsonObject.containsKey("mature_mirnas_file")) {
                Path valid_mirnas_file_Path = Paths.get(jsonObject.getString("mature_mirnas_file"));
                if (Files.exists(valid_mirnas_file_Path)) {
                    mirbaseFile = valid_mirnas_file_Path.toFile();
                    System.out.println("Mature miRNAs FASTA file (from miRBase) is loaded.");
                } else {
                    throw new IOException("Mature miRNAs FASTA file (from miRBase): " + valid_mirnas_file_Path + " not found.");
                }
            }else{
                throw new IOException("Mature miRNAs FASTA file (from miRBase) is required.");
            }
        }else if (jsonObject.containsKey("mature_mirnas_file")) {
            Path valid_mirnas_file_Path = Paths.get(jsonObject.getString("mature_mirnas_file"));
            if (Files.exists(valid_mirnas_file_Path)) {
                mirbaseFile = valid_mirnas_file_Path.toFile();
                System.out.println("Mature miRNAs FASTA file (from miRBase) is loaded.");
            } else {
                throw new IOException("Mature miRNAs FASTA file (from miRBase): " + valid_mirnas_file_Path + " not found.");
            }
        }
        
        if (jsonObject.containsKey("output_dir")) {
            Path output_dir_path = Paths.get(jsonObject.getString("output_dir"));
            if (!Files.exists(output_dir_path)) {
                output_dir_path.toFile().mkdirs();
            }
            outputDir = output_dir_path.toFile();
        } else {
            //output current working directory message
            Path currentRelativePath = Paths.get("");
            String s = currentRelativePath.toAbsolutePath().toString();
            System.out.println("Output directory not set. Using current working directory: " + s + File.separator + "PAREnet");
            outputDir = new File(s + File.separator + "PAREnet");
            if (outputDir.exists() && outputDir.listFiles().length > 0) {
                System.out.println("Default generated directory: " + s + File.separator + "PAREnet" + " exists and is not empty. Exiting...");
                System.exit(1);
            }
            outputDir.mkdir();
            
            if (new File(outputDir.getAbsolutePath() + File.separator + "PAREnet_CysIn.csv").exists()) {
                System.out.println("A PAREnet output file already exists in this directory.\nExiting...");
                System.exit(1);
            }

        }
                */
    }
}
