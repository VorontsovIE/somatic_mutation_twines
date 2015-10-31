require 'WingenderTFClass'
require_relative 'ensembl_mappings'

HgncByEnsembl = load_hgnc_by_ensembl('source_data/EnsemblToHGNC_GRCh37.p13.tsv')

Twine = Struct.new(:chromosome, :interval, :sites, :ensembl_gene_ids) do
  def hgnc_names
    @hgnc_names ||= ensembl_gene_ids.flat_map{|ensembl_gene_id| HgncByEnsembl[ensembl_gene_id] }.uniq.compact
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
  
  def mutations
    sites.map(&:variant_id).uniq
  end

  def samples
    sites.map(&:sample_name).uniq
  end

  def num_mutations
    mutations.size
  end

  def num_samples
    samples.size
  end

  def length
    interval.to_interval_set.total_length
  end
  
  def to_s
    '<' + \
    [
      "chr#{chromosome}:#{interval.integer_points.first}-#{interval.integer_points.last}", 
      "#{num_samples} samples", 
      "#{num_mutations} mutations: " '{' + mutations.join(', ') + '}', 
      # '{' + involved_families.join('; ') + '}',
    ].join('; ') + \
    '>'
  end
  def inspect; to_s; end
end
