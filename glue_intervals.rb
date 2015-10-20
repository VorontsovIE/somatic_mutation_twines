$:.unshift File.absolute_path('./lib', __dir__)
require_relative 'configuration'
require 'snv_info'
require 'perfectosape/SNPScanResults'
require 'interval_notation'
require 'WingenderTFClass'

EnsemblMappings = Struct.new(:ensg, :enst, :hgnc_symbol, :description, :uniprot_id, :uniprot_ac) do
  def self.from_string(str)
    ensg, enst, hgnc_symbol, description, uniprot_id, uniprot_ac = str.chomp.split("\t", 6)
    self.new(ensg, enst, hgnc_symbol, description, uniprot_id, uniprot_ac)
  end

  def self.each_in_file(filename, &block)
    File.readlines(filename).map{|line|
      self.from_string(line)
    }.each(&block)
  end
end

def load_hgnc_by_ensembl(mapping_filename)
  EnsemblMappings.each_in_file(mapping_filename).group_by(&:ensg).map{|ensg, infos|
    [ensg, infos.map(&:hgnc_symbol).uniq.compact]
  }.to_h
end

def load_uniprot_id_by_ensembl(mapping_filename)
  EnsemblMappings.each_in_file(mapping_filename).group_by(&:ensg).map{|ensg, infos|
    [ensg, infos.map(&:uniprot_id).uniq.compact]
  }.to_h
end

HgncByEnsembl = load_hgnc_by_ensembl('source_data/EnsemblToHGNC_GRCh37.p13.tsv')
UniprotByEnsembl = load_uniprot_id_by_ensembl('source_data/EnsemblToHGNC_GRCh37.p13.tsv')



Twine = Struct.new(:chromosome, :interval, :mutations, :sites, :ensembl_gene_ids) do
  def hgnc_names
    @hgnc_names ||= ensembl_gene_ids.map{|ensembl_gene_id| HgncByEnsembl[ensembl_gene_id] }.uniq.compact
  end
  def involved_motifs
    sites.map(&:motif_name).uniq
  end
  def involved_families
    involved_motifs.flat_map{|motif|
      uniprot_id = motif.to_s.split('.').first
      WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[3].subfamilies_by_uniprot_id(uniprot_id)
    }.uniq
  end
  def to_s
    '<' + ["chr#{chromosome}", interval, '{' + mutations.map(&:position).join(',') + '}', ensembl_gene_ids, hgnc_names].join('; ') + '>'
  end
  def inspect; to_s; end
end

chromosomes = (1..22).map(&:to_s).map(&:to_sym) + [:X, :Y]

cancer_type = ARGV[0]
dhs_accessible_filename = ARGV[1]
sites_filename = "results/sites/#{cancer_type}.txt"
mutations_filename = "source_data/SNVs/#{cancer_type}.txt"



# markup_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty)
# puts main_chromosomes(GENOME_READER).map{|chr| ensembl_markup_by_chromosome[chr] }.map(&:total_length).inject(&:+)

mutations_by_chromosome = SNVInfo
  .each_in_file(mutations_filename).to_a
  .group_by(&:chromosome)

sites_by_chromosome = PerfectosAPE::Result
  .each_in_file(sites_filename).to_a
  .group_by(&:chromosome).map{|chromosome, sites|
    [chromosome, sites.group_by(&:variant_id)]
  }.to_h


twines_by_chromosome = PerfectosAPE::Result
  .each_in_file(sites_filename).to_a
  .group_by(&:chromosome).map{|chromosome, sites|
    intervals = sites.map(&:strongest_site_interval)
    [chromosome, IntervalNotation::Operations.union(intervals)]
  }.to_h
twines_by_chromosome.default = IntervalNotation::Empty


genes_by_chromosome = load_genes_by_chromosome('source_data/ensembl_markup.tsv', length_5_prime: 5000, length_3_prime: 0)

twines = chromosomes.flat_map do |chromosome|
  mutations = mutations_by_chromosome[chromosome]
  sites = sites_by_chromosome[chromosome]
  gene_intervals_by_name = genes_by_chromosome[chromosome]
  twines_by_chromosome[chromosome].intervals.map(&:to_interval_set).map{|twine_interval|
    selected_mutations = mutations.select{|mutation| twine_interval.include_position?(mutation.position) }
    ensembl_gene_ids = gene_intervals_by_name.select{|ensembl_gene_id, interval|
      interval.intersect?(twine_interval)
    }.keys
    ensembl_gene_ids = []
    selected_sites = selected_mutations.flat_map{|mutation| sites[mutation.variant_id] }
    Twine.new(chromosome, twine_interval, selected_mutations, selected_sites, ensembl_gene_ids)
  }
end
puts twines.select{|twine| twine.mutations.size > 1 }.first(10).map{|twine|
  [twine, twine.involved_motifs, twine.involved_families].inspect
}



# chromosomes.reject{|chr|
#   num_mutations_by_chromosome[chr] == 0 || twines_by_chromosome[chr].empty?
# }.sort_by{|chr|
#   # twines_by_chromosome[chr].total_length.to_f / num_mutations_by_chromosome[chr]
#   num_mutations_by_chromosome[chr].to_f / twines_by_chromosome[chr].num_connected_components
# }.reverse_each do |chr|
#   intervals = twines_by_chromosome[chr]
#   num_mutations = num_mutations_by_chromosome[chr]

#   aveComponentLength = intervals.total_length.to_f / intervals.num_connected_components
#   mutDivLen = num_mutations.to_f / intervals.total_length
#   mutDivComp = num_mutations.to_f / intervals.num_connected_components
#   puts "chr#{chr}\tcomp: #{intervals.num_connected_components};\tlength: #{intervals.total_length}\tmutation: #{num_mutations}\taveComponentLength: #{aveComponentLength.round(1)}\tmutDivLen: #{mutDivLen.round(2)}\tmutDivComp: #{mutDivComp.round(2)}"
# end