$:.unshift File.absolute_path('./lib', __dir__)
require_relative 'configuration'
require 'set'
require 'snv_info'
require 'perfectosape/SNPScanResults'
require 'interval_notation'
require 'yaml'
require 'ensembl_mappings'
require 'twine'

# length of regulatory region of a gene
def read_gene_length(filename)
  File.readlines(filename).map(&:chomp).map{|l|
    l.split("\t")
  }.map{|hgnc, length|
    [hgnc, length.to_i]
  }.to_h
end


raise 'Specify twines filename'  unless twines_filename = ARGV[0]
length_by_gene = read_gene_length('regulatory_gene_length.tsv')
twines_by_gene = YAML.load_file(twines_filename)


HgncByEnsembl = load_hgnc_by_ensembl('source_data/EnsemblToHGNC_GRCh37.p13.tsv')
protein_coding_genes = File.readlines('source_data/allhgnc_protein_coding.txt').map(&:chomp).to_set
hgncs_all = HgncByEnsembl.values.flatten.uniq.to_set # & protein_coding_genes
# oncogenes = File.readlines('source_data/сancer_genes_lists/allonco_20130923_reduced.csv')
#   .map{|l|
#     l.chomp.split
#   }.map{|gene, num_evidences|
#     [gene, num_evidences.to_i]
#   }
# reliable_oncogenes = oncogenes
#   .select{|gene, num_evidences|
#     num_evidences >= 1
#   }.map{|gene, num_evidences|
#     gene
#   }.to_set
reliable_oncogenes = File.readlines('source_data/сancer_genes_lists/onco_gathered.txt').map(&:chomp).to_set

known_oncogenes = reliable_oncogenes & hgncs_all



sorted_genes_with_twines = twines_by_gene.sort_by{|gene, twines|
  twines.map(&:num_mutations).max # most useful (place 1)
  # twines.map(&:num_samples).max # useful (place 2)
}.reverse

puts '#' + [
  'Gene',
  'Known oncogene',
  'Num twines',
  'Regulatory length',
  'Max mutations in twine',
  'Max unique mutations in twine',
  'Max samples in twine',
  'Max twine length',
  'Max density of mutations in twine',
  'Total num of mutations',
  'Total length of twines',
].join("\t")

sorted_genes_with_twines.each{|gene, twines|
  puts [
    gene,
    known_oncogenes.include?(gene) ? 1 : 0,
    twines.size,
    length_by_gene[gene],
    twines.map(&:num_mutations).max,
    twines.map(&:num_unique_mutations).max,
    twines.map(&:num_samples).max,
    twines.map(&:length).max,
    twines.map(&:mutation_density).max,
    twines.map(&:num_mutations).inject(0, &:+),
    twines.map(&:length).inject(0, &:+),
  ].join("\t")
}
