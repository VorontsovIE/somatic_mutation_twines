query_xml = <<-EOS
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
            
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "hgnc_symbol" />
        <Attribute name = "description" />
        <Attribute name = "uniprot_swissprot" />
        <Attribute name = "uniprot_swissprot_accession" />
    </Dataset>
</Query>
EOS

desc 'Download and prepare Ensembl to HGNC/to Uniprot mapping (Ensembl, GRCh37.p13)'
task :download_mappings => 'source_data/EnsemblToHGNC_GRCh37.p13.tsv'
file 'source_data/EnsemblToHGNC_GRCh37.p13.tsv' do
  sh 'wget', '-O', 'source_data/EnsemblToHGNC_GRCh37.p13.tsv', "http://feb2014.archive.ensembl.org/biomart/martservice?query=#{query_xml}"
end
