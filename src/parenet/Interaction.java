/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package parenet;

import java.util.Objects;


/**
 *
 * @author twh14ura
 */
class Interaction {

    private String srna_node_color = "#8ED1E3";	//light blue    #EDC39A";    // light orange     
    private final String srna_seq;
    private String sRNA_node_label;
    private String mRNA_node_label;
    private String miRid;
   
    private int category;
    private final int cleavePos;
    private boolean isValidatedInteraction;
    private boolean isKnownMirna;
    private final boolean isSourceEdge;
    private final Transcript mRNA;
                
    public void setMirnaAnnotation(String miRid) {
        this.miRid = miRid;
        this.sRNA_node_label = miRid;
        this.isKnownMirna = true;
        srna_node_color = "#EDC39A";    // light orange
    }

    public void setValidated() {
        this.isValidatedInteraction = true;
    }
    
    public void setCategory(int category) {
        this.category = category;
    }

    public int getCategory() {
        return category;
    }
    
    public Interaction(String srna_seq, int category, int cleavePos, Transcript mRNA) {
        this.srna_seq = srna_seq;
        this.miRid = srna_seq;
        this.sRNA_node_label = "";
        this.mRNA_node_label = mRNA.getGene_name().split(";")[0].trim();
        this.category = category;
        this.cleavePos = cleavePos;
        this.isKnownMirna = false;
        this.isValidatedInteraction = false;
        this.isSourceEdge = false;
        this.mRNA = mRNA;
    }
    
    public Interaction(String srna_seq, String miRid, boolean isKnownMirna, Transcript mRNA){
        this.srna_seq = srna_seq;
        this.miRid = miRid;
        this.sRNA_node_label = "";
        this.mRNA_node_label = mRNA.getGene_name().split(";")[0].trim();
        this.category = -1;
        this.cleavePos = -1;
        this.isKnownMirna = isKnownMirna;
        this.isValidatedInteraction = false;
        this.isSourceEdge = true;
        this.mRNA = mRNA;
    }

    public String getMiRid() {
        return miRid;
    }

    public boolean isIsKnownMirna() {
        return isKnownMirna;
    }

    public Transcript getmRNA() {
        return mRNA;
    }
    
    public String getGene_ID() {
        return mRNA.getGene_id();
    }
    
    public String getGene_name(){
        return mRNA.getGene_name();
    }

    @Override
    public String toString() {
        return miRid + "," + srna_seq + "," + isKnownMirna + "," + sRNA_node_label + "," + srna_node_color + "," + mRNA.getGene_id() + ",\"" + mRNA.getGene_description() + "\"," + mRNA_node_label + "," + mRNA.getGene_node_color()  + "," +  mRNA.isSource() + "," +  isSourceEdge + "," + isValidatedInteraction + "," + category + "," + cleavePos;
    }
    
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 47 * hash + Objects.hashCode(this.mRNA.getGene_id());
        hash = 47 * hash + (this.isSourceEdge ? 1 : 0);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Interaction other = (Interaction) obj;
        if (!Objects.equals(this.mRNA.getGene_id(), other.mRNA.getGene_id())) {
            return false;
        }
        if (this.isSourceEdge != other.isSourceEdge) {
            return false;
        }
        return true;
    }

}
