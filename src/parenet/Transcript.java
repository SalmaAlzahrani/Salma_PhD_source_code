/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package parenet;

/**
 *
 * @author twh14ura
 */
class Transcript {

    private final String gene_id;
    private final String gene_description;
    private String gene_name;
    private String gene_node_color = "#CCCCFF"; // purple
    private boolean isSource;

    public Transcript(String gene_id, String gene_name, String gene_description) {
        this.gene_id = gene_id;
        this.gene_name = gene_name;
        this.gene_description = gene_description;
        this.isSource = false;
    }

    public void setSource(){
        this.isSource = true;
    }
    
    // change color of node if it's annotated gene
    public void setGeneName(String geneName){
        this.gene_name = geneName;
        this.gene_node_color = "#7FCDBB"; // green
    }
    
    // do sRNA originate from this gene
    public boolean isSource(){
        return isSource;
    }
   
    public String getGene_id() {
        return gene_id;
    }

    public String getGene_description() {
        return gene_description;
    }

    public String getGene_name() {
        return gene_name;
    }

    public String getGene_node_color() {
        return gene_node_color;
    }

//    public boolean isIsSource() {
//        return isSource;
//    }
        
}
