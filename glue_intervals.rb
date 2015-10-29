$:.unshift File.absolute_path('./lib', __dir__)
require_relative 'configuration'
require 'set'
require 'snv_info'
require 'perfectosape/SNPScanResults'
require 'interval_notation'
require 'yaml'
require 'ensembl_mappings'
require 'twine'

# def load_uniprot_id_by_ensembl(mapping_filename)
#   EnsemblMappings.each_in_file(mapping_filename).group_by(&:ensg).map{|ensg, infos|
#     [ensg, infos.map(&:uniprot_id).uniq.compact]
#   }.to_h
# end

# UniprotByEnsembl = load_uniprot_id_by_ensembl('source_data/EnsemblToHGNC_GRCh37.p13.tsv')



def gene_twines(twines)
  result = Hash.new{|h, k| h[k] = [] }
  twines.each{|twine|
    twine.hgnc_names.reject(&:nil?).reject(&:empty?).each{|hgnc|
      result[hgnc] << twine
    }
  }
  result
end

def make_segmentations_by_chromosome(genes_by_chromosome)
  genes_by_chromosome.map{|chromosome, gene_by_name|
    gene_intervals_tagged = gene_by_name.map{|name, interval| [interval, name.to_s] }
    segmentation = IntervalNotation::SweepLine.make_tagging(gene_intervals_tagged)
    [chromosome, segmentation]
  }.to_h
end

def load_sites_by_chromosome(sites_filename)
  result = PerfectosAPE::Result
    .each_in_file(sites_filename).to_a
    .group_by(&:chromosome)
  result.default_proc = ->(hsh, chromosome){ [] }
  result
end

def make_site_union_by_chromosome(sites_by_chromosome)
  result = sites_by_chromosome.map{|chromosome, sites|
    intervals = sites.map(&:strongest_site_interval)
    [chromosome, IntervalNotation::Operations.union(intervals)]
  }.to_h
  result.default = IntervalNotation::Empty
  result
end


chromosomes = (1..22).map(&:to_s).map(&:to_sym) + [:X, :Y]

raise 'Specify cnacer type'  unless cancer_type = ARGV[0]
sites_filename = "results/sites/#{cancer_type}.txt"
# mutations_filename = "source_data/SNVs/#{cancer_type}.txt"


tm = Time.now

genes_by_chromosome = load_genes_by_chromosome('source_data/ensembl_markup.tsv', length_5_prime: 5000, length_3_prime: 0)

$stderr.puts "load gene markup #{Time.now - tm}"; tm = Time.now

segmentations_by_chromosome = make_segmentations_by_chromosome(genes_by_chromosome)

$stderr.puts "make gene segmentation #{Time.now - tm}"; tm = Time.now

sites_by_chromosome = load_sites_by_chromosome(sites_filename)

$stderr.puts "load sites #{Time.now - tm}"; tm = Time.now

site_union_by_chromosome = make_site_union_by_chromosome(sites_by_chromosome)

$stderr.puts "sites union #{Time.now - tm}"; tm = Time.now

twines = chromosomes.flat_map do |chromosome|
  gene_segmentation = segmentations_by_chromosome[chromosome]
  gene_intervals_by_name = genes_by_chromosome[chromosome]

  sites_by_twine = Hash.new{|h,k| h[k] = [] }
  genes_by_twine = Hash.new{|h,k| h[k] = Set.new }

  sites = sites_by_chromosome[chromosome]

  site_intervals = site_union_by_chromosome[chromosome]

  sites.each do |site|
    twine_interval = site_intervals.interval_covering_point(site.mutation_genomic_position)
    sites_by_twine[twine_interval] << site

    genes_on_mutation = gene_segmentation.segment_covering_point(site.mutation_genomic_position).state
    genes_by_twine[twine_interval] += genes_on_mutation
  end

  site_intervals.intervals.map{|twine_interval|
    Twine.new(chromosome, twine_interval, sites_by_twine[twine_interval], genes_by_twine[twine_interval])
  }.select{|twine|
    twine.mutations.size > 1 && twine.samples.size > 1
  }
end

$stderr.puts "aggregate twines #{Time.now - tm}"; tm = Time.now


twines_ensembl = twines.select{|twine| !twine.ensembl_gene_ids.empty? }
twines_hgnc = twines.select{|twine| !twine.hgnc_names.empty? }

twines_by_gene = gene_twines(twines_hgnc)

puts twines_by_gene.to_yaml
