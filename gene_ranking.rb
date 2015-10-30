$:.unshift File.absolute_path('./lib', __dir__)
require_relative 'configuration'
require 'set'
require 'snv_info'
require 'perfectosape/SNPScanResults'
require 'interval_notation'
require 'yaml'
require 'ensembl_mappings'
require 'twine'

raise 'Specify cancer type'  unless cancer_type = ARGV[0]

# length of regulatory region of a gene
def read_gene_length(filename)
  File.readlines(filename).map(&:chomp).map{|l|
    l.split("\t")
  }.map{|hgnc, length|
    [hgnc, length.to_i]
  }.to_h
end

length_by_gene = read_gene_length('regulatory_gene_length.tsv')
# full_length_by_gene = read_gene_length('gene_length.tsv')

twines_by_gene = YAML.load_file("results/twines/#{cancer_type}.yaml")

# puts twines_by_gene.size

twines_by_gene.sort_by{|gene, twines|
  # twines.flat_map(&:samples).uniq.size.to_f / length_by_gene[gene]
  twines.size.to_f / length_by_gene[gene]
  # twines.map(&:num_samples).inject(0.0, &:+) / length_by_gene[gene]
  # twines.flat_map(&:mutations).size.to_f / twines.size
  # twines.flat_map(&:mutations).size.to_f / twines.map(&:length).inject(0.0, &:+)
  # twines.map(&:length).inject(0.0, &:+) / length_by_gene[gene]
  # twines.flat_map{|twine| twine.mutations.size.to_f / twine.length }.max
  # twines.flat_map(&:mutations).size / twines.flat_map(&:length).inject(&:+).to_f
  # twines.map(&:length).max
  # twines.map(&:num_samples).max # most useful
  # twines.map(&:num_mutations).max # biased for kataegis
}.each{|gene, twines|
  total_samples = twines.flat_map(&:samples).uniq.size
  total_samples_distinct = twines.flat_map(&:samples).size
  total_mutations = twines.flat_map(&:mutations).uniq.size
  puts "#{gene} (#{length_by_gene[gene]} bp) -- #{twines.size} twines -- #{total_mutations} mutations -- #{ total_samples } (#{total_samples_distinct}) samples (distinct by twine)"
  puts twines
  puts "num bp per twine #{(length_by_gene[gene].to_f / twines.size).round(1) }"
  puts '---------------------------------'
}

all_twines = twines_by_gene.flat_map{|gene, twines| twines }
all_twines_ensembl = all_twines.select{|twine| !twine.ensembl_gene_ids.empty? }
all_twines_hgnc = all_twines.select{|twine| !twine.hgnc_names.empty? }


puts
puts({
  twines_by_gene: twines_by_gene.size,
  multitwines_by_gene: twines_by_gene.count{|gene, twines| twines.size > 1 },
  twines: all_twines.size,
  twines_ensembl: all_twines_ensembl.size,
  twines_hgnc: all_twines_hgnc.size,
})
