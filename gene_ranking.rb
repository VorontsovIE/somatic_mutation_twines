$:.unshift File.absolute_path('./lib', __dir__)
require_relative 'configuration'
require 'set'
require 'snv_info'
require 'perfectosape/SNPScanResults'
require 'interval_notation'
require 'yaml'
require 'ensembl_mappings'
require 'twine'

raise 'Specify cnacer type'  unless cancer_type = ARGV[0]

twines_by_gene = YAML.load_file("results/twines/#{cancer_type}.yaml")

# puts twines_by_gene.size

twines_by_gene.select{|gene, twines|
  twines.size > 1
}.sort_by{|gene, twines|
  -twines.size
}.each{|gene, twines|
  $stderr.puts "#{gene} -- #{twines.size} twines"
  $stderr.puts twines
  $stderr.puts '---------------------------------'
}

all_twines = twines_by_gene.flat_map{|gene, twines| twines }
all_twines_ensembl = all_twines.select{|twine| !twine.ensembl_gene_ids.empty? }
all_twines_hgnc = all_twines.select{|twine| !twine.hgnc_names.empty? }


$stderr.puts
$stderr.puts({
  twines_by_gene: twines_by_gene.size,
  multitwines_by_gene: twines_by_gene.count{|gene, twines| twines.size > 1 },
  twines: all_twines.size,
  twines_ensembl: all_twines_ensembl.size,
  twines_hgnc: all_twines_hgnc.size,
})
